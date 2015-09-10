// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/connection/BridgeChains.cc
/// @brief The BridgeChains mover, for connecting chains with fixed relative orientation
/// @details
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/connection/BridgeChains.hh>
#include <protocols/denovo_design/connection/BridgeChainsCreator.hh>

//Project Headers

//Core Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constants.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/sequence/ABEGOManager.hh>

//Protocol Headers
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/DsspMover.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers
#include <math.h>

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.denovo_design.BridgeChains" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace connection {
////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace ObjexxFCL;

std::string
BridgeChainsCreator::keyname() const
{
	return BridgeChainsCreator::mover_name();
}

protocols::moves::MoverOP
BridgeChainsCreator::create_mover() const
{
	return protocols::moves::MoverOP( new BridgeChains() );
}

std::string
BridgeChainsCreator::mover_name()
{
	return "BridgeChains";
}

///  ---------------------------------------------------------------------------------
///  BridgeChains main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
BridgeChains::BridgeChains() :
	protocols::moves::Mover( "BridgeChains" ),
	motif_( "" ),
	chain1_( 1 ),
	chain2_( 2 ),
	overlap_( 3 ),
	scorefxn_( /* NULL */ ),
	vlb_( protocols::forge::components::VarLengthBuildOP( new protocols::forge::components::VarLengthBuild() ) ),
	cached_ss_( "" ),
	cached_aa_( "" ),
	cached_start_( 0 ),
	cached_end_( 0 )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
BridgeChains::~BridgeChains()
{}


/// Return a copy of ourselves
protocols::moves::MoverOP
BridgeChains::clone() const
{
	return protocols::moves::MoverOP( new BridgeChains(*this) );
}

std::string
BridgeChains::get_name() const {
	return BridgeChainsCreator::mover_name();
}

void
BridgeChains::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	motif_ = tag->getOption< std::string >( "motif", motif_ );
	chain1_ = tag->getOption< core::Size >( "chain1", chain1_ );
	chain2_ = tag->getOption< core::Size >( "chain2", chain2_ );
	overlap_ = tag->getOption< core::Size >( "overlap", overlap_ );
	std::string const sfxn ( tag->getOption<std::string>( "scorefxn", "" ) );
	if( sfxn != "" ) {
		scorefxn_ = data.get_ptr<core::scoring::ScoreFunction>( "scorefxns", sfxn );
		TR << "score function " << sfxn << " is used. " << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void BridgeChains::apply( core::pose::Pose & pose )
{
	runtime_assert( scorefxn_ != 0 );
	// DSSP the heck out of this pose
	protocols::moves::DsspMover dssp;
	dssp.apply( pose );

	// check jump input to make sure it's valid
	runtime_assert( chain1_ <= pose.conformation().num_chains() );
	runtime_assert( chain2_ <= pose.conformation().num_chains() );
	// left is the c-terminal residue of jump1
	core::Size const start( pose.conformation().chain_end( chain1_ ) );
	// right is the n-terminal residue of jump2
	core::Size const end( pose.conformation().chain_begin( chain2_ ) );
	// use overlap to expand rebuilt area
	core::Size left( start - overlap_ );
	core::Size right( end + overlap_ );
	if ( left < 1 ) {
		left = 1;
	} else if ( left > pose.total_residue() ) {
		left = pose.total_residue();
	}
	if ( right < 1 ) {
		right = 1;
	} else if ( right > pose.total_residue() ) {
		right = pose.total_residue();
	}

	utility::vector1< std::string > const & motifs( utility::string_split( motif_, '-' ) );
	core::Size insert_length( 0 );
	std::string ss( "" );
	std::string aa( "" );
	utility::vector1< std::string > abego_insert;
	core::scoring::constraints::ConstraintCOPs constraint_set;
	// residues to be connected are obliterated by the vlb, so we need to add them to the ss/abego
	for ( core::Size i=left; i<=start; ++i ) {
		aa += pose.residue( i ).name1();
		abego_insert.push_back( "X" );
		ss += pose.secstruct( i );
		constraint_set.push_back( create_coordinate_cst( pose, i ) );
	}
	for ( core::Size i=1; i<=motifs.size(); ++i ) {
		// here, we can accept "3LX" or "3:LX"
		std::string motif( "" );
		for ( core::Size j=0; j<motifs[i].size(); ++j ) {
			if ( motifs[i][j] != ':' ) {
				motif += motifs[i][j];
			}
		}
		std::string const ss_type( motif.substr( motif.size()-2, 1 ) );
		std::string const abego_type( motif.substr( motif.size()-1, 1 ) );
		assert( ss_type.size() );
		assert( abego_type.size() );
		if ( ss_type[0] != 'E' && ss_type[0] != 'H' && ss_type[0] != 'L' ) {
			TR.Error << "Invalid SS type in motif " << motif << std::endl;
			utility_exit();
		}
		if ( abego_type[0] > 'Z' || abego_type[0] < 'A' ) {
			TR.Error << "Invalid abego type in motif " << motif << std::endl;
			utility_exit();
		}
		core::Size const len( (core::Size)utility::string2int( motif.substr( 0, motif.size()-2 ) ) );
		TR << "motif" << i << " = " << motif << " " << ss_type << " " << abego_type << " " << len << std::endl;
		insert_length += len;
		for ( core::Size j=1; j<=len; ++j ) {
			ss += ss_type;
			aa += "V";
			abego_insert.push_back( abego_type );
		}
	}
	for ( core::Size i=end; i<=right; ++i ) {
		// residues on the right must also be added
		aa += pose.residue( i ).name1();
		abego_insert.push_back( "X" );
		ss += pose.secstruct( i );
		constraint_set.push_back( create_coordinate_cst( pose, i ) );
	}
	TR << "Original size is " << pose.total_residue() << "; Going to build the secondary structure: " << ss << " between residues " << left << " and " << right << " to connect the jump starting at residue " << start << " and ending at residue " << end << std::endl;

	utility::vector1< std::string > abego;
	runtime_assert( left < right );
	for ( core::Size i=1; i<left; ++i ) {
		abego.push_back( "X" );
	}
	for ( core::Size j=1; j<=abego_insert.size(); ++j ) {
		abego.push_back( abego_insert[ j ] );
	}

	// setup complete, now start working with the pose

	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;

	runtime_assert( vlb_ != 0 );
	if ( ( ss != cached_ss_ ) || ( aa != cached_aa_ ) || ( cached_start_ != left ) || ( cached_end_ != right ) ) {
		// the build manager
		protocols::forge::build::BuildManager manager;

		// add the instruction to rebuild this segment
		manager.add( protocols::forge::build::BuildInstructionOP( new protocols::forge::build::SegmentRebuild( protocols::forge::build::Interval( left, right ), ss, aa ) ) );

		//clear fragment cache and set build manager
		vlb_->manager( manager );

		// set cached values
		cached_ss_ = ss;
		cached_aa_ = aa;
		cached_start_ = left;
		cached_end_ = right;
	}
	vlb_->scorefunction( scorefxn_->clone() );
	vlb_->vall_memory_usage( protocols::forge::components::VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->loop_mover_str( "RemodelLoopMover" );
	vlb_->set_abego( abego );

	// add coordinate constraint generator
	//protocols::forge::remodel::RemodelConstraintGeneratorOP cst_gen(
	//		new protocols::forge::remodel::RemodelConstraintGenerator() );
	// add coordinate constraints for all overlapping residues
	//cst_gen->add_constraints( constraint_set );
	//vlb_->clear_rcgs();
	//vlb_->add_rcg( cst_gen );

	bool closed( false );
	core::pose::Pose test_pose( modified_archive_pose );
	vlb_->apply( test_pose );
	if ( vlb_->get_last_move_status() == protocols::moves::MS_SUCCESS ) {
		closed = true;
	}
	if ( closed ) {
		TR << "CLOSED THE LOOP" << std::endl;
		set_last_move_status( protocols::moves::MS_SUCCESS );
		pose = test_pose;
		TR << "connected abego = " << core::sequence::get_abego( pose, 2 ) << std::endl;
		// check SS against wanted SS
		dssp.apply( pose );
		std::string pose_ss( pose.secstruct() );
		bool match( true );
		for ( core::Size i=1; i<=ss.size(); ++i ) {
			if ( pose_ss[i+left-2] != ss[i-1] ) {
				match = false;
				TR << "Connection in pose doesn't match the desired secondary structure. Wanted= " << ss << "; Actual= " << pose_ss << std::endl;
				break;
			}
		}
		if ( !match ) {
			test_pose.dump_pdb( "ss_failed.pdb" );
			set_last_move_status( protocols::moves::FAIL_RETRY );
		}
	} else {
		test_pose.dump_pdb( "failed.pdb" );
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}
}

/// @brief creates a Ca coordinate constraint for residue resi
core::scoring::constraints::ConstraintOP
BridgeChains::create_coordinate_cst( core::pose::Pose const & pose,
		core::Size const resi ) const
{
	core::Size atom( pose.residue_type(resi).nbr_atom() );
	if( pose.residue_type(resi).has("CA") ) {
		atom = pose.residue_type(resi).atom_index("CA");
	}

	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc(0.0, 0.5) );
	return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint(
			core::id::AtomID(atom,resi),
			core::id::AtomID(pose.residue(1).nbr_atom(),1),
			pose.residue(resi).xyz(atom),
			fx ) );
}

} // namespace connection
} // namespace denovo_design
} // namespace protocols