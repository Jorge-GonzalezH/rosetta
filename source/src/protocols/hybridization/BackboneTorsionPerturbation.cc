// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Yifan Song

#include <protocols/hybridization/BackboneTorsionPerturbation.hh>
#include <protocols/hybridization/BackboneTorsionPerturbationCreator.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/import_pose/import_pose.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/hybridization/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/USOGFunc.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <ObjexxFCL/format.hh>

// task operation
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/random/random.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/datacache/DataMap.hh> // AUTO IWYU For DataMap

static basic::Tracer TR( "protocols.hybridization.BackboneTorsionPerturbation" );

namespace protocols {
namespace hybridization {

using namespace core;

BackboneTorsionPerturbation::BackboneTorsionPerturbation() {
	init();
}

BackboneTorsionPerturbation::~BackboneTorsionPerturbation() = default;

void BackboneTorsionPerturbation::init() {
	increase_cycles_ = 1.0;
	temperature_ = 5.0;
	recover_low_ = false;
	local_ = 0;
	dump_snapshots_ = false;
	snapshot_interval_ = 10;
}

void add_constraints(core::pose::Pose & pose, core::Size rsd1, core::Size rsd2) {
	using namespace core::scoring::func;
	TR << "constrain " << rsd1 << " " << rsd2 << std::endl;
	for ( core::Size ires=1; ires<rsd1-1; ++ires ) {
		if ( !pose.residue_type(ires).is_protein() ) continue;
		for ( core::Size jres=rsd2+2; jres<=pose.size(); ++jres ) {
			if ( !pose.residue_type(jres).is_protein() ) continue;

			core::Real COORDDEV = 1.0;
			core::Real dist = pose.residue(ires).xyz(2).distance( pose.residue(jres).xyz(2) );
			FuncOP fx( new ScalarWeightedFunc( 1.0, utility::pointer::make_shared< USOGFunc >( dist, COORDDEV ) ) );
			pose.add_constraint(
				utility::pointer::make_shared< core::scoring::constraints::AtomPairConstraint >( core::id::AtomID(2,ires), core::id::AtomID(2,jres), fx )
			);
		}
	}
}

void optimize(core::pose::Pose & pose, core::Size rsd1, core::Size rsd2, core::scoring::ScoreFunctionOP scorefxn, core::Size ncycles, core::Real max_delta_torsion) {
	protocols::moves::MonteCarloOP mc;
	mc = utility::pointer::make_shared< protocols::moves::MonteCarlo >( pose, *scorefxn, 2 );

	for ( core::Size icycle=1; icycle<=ncycles; ++ icycle ) {
		core::Size ires = numeric::random::rg().random_range(rsd1, rsd2);
		core::Real delta = (2.* numeric::random::rg().uniform() - 1.) * max_delta_torsion;
		if ( numeric::random::rg().uniform() < 0.5 ) {
			pose.set_phi(ires, pose.phi(ires) + delta);
		} else {
			pose.set_psi(ires, pose.psi(ires) + delta);
		}
		//core::Real score=(*scorefxn)(pose);
		mc->boltzmann(pose, "Optimization");
	}
	mc->show_scores();
	mc->show_counters();
	mc->recover_low(pose);
}

void pick_pivots(core::pose::Pose & pose, core::Size & rsd1, core::Size & torsion1, core::Size & rsd2, core::Size & torsion2, core::Size variance=5) {
	// pick two torsion as pivots
	core::Size n_torsion(0);
	for ( core::Size ires=1; ires <= pose.size(); ++ires ) {
		n_torsion += 2; // for now
	}

	int perturbed_torsion(0), coupled_torsion(0);
	while ( coupled_torsion <= 0 || coupled_torsion > (int) n_torsion ) {
		perturbed_torsion = numeric::random::rg().random_range(1, n_torsion);
		rsd1 = (perturbed_torsion + 1)/2;
		torsion1 = perturbed_torsion - (rsd1-1)*2;

		coupled_torsion = perturbed_torsion;
		while ( coupled_torsion == (int) perturbed_torsion ) {
			coupled_torsion = floor( perturbed_torsion + variance*numeric::random::rg().gaussian() + 0.5 );
		}
		rsd2 = (coupled_torsion + 1)/2;
		torsion2 = coupled_torsion - (rsd2-1)*2;
	}

	TR << perturbed_torsion << " " << coupled_torsion << std::endl;
}

void BackboneTorsionPerturbation::perturb(core::pose::Pose & pose,
	core::Real max_delta_torsion) {

	core::Real delta = (2.* numeric::random::rg().uniform() - 1.) * max_delta_torsion;
	core::Size rsd1, rsd2, torsion1, torsion2;
	pick_pivots(pose, rsd1, torsion1, rsd2, torsion2);

	if ( rsd1 < rsd2 ) {
		add_constraints(pose, rsd1, rsd2);
	} else {
		add_constraints(pose, rsd2, rsd1);
	}

	if ( torsion1 == 1 ) {
		pose.set_phi(rsd1, pose.phi(rsd1) + delta);
	} else {
		pose.set_psi(rsd1, pose.psi(rsd1) + delta);
	}
	if ( torsion2 == 1 ) {
		pose.set_phi(rsd2, pose.phi(rsd2) - delta);
	} else {
		pose.set_psi(rsd2, pose.psi(rsd2) - delta);
	}

	core::scoring::ScoreFunctionOP cst_fxn = scorefxn_->clone();
	cst_fxn->reset();
	cst_fxn->set_weight( core::scoring::atom_pair_constraint, 1. );
	optimize(pose, rsd1, rsd2, cst_fxn, 100, 5.);

	pose.remove_constraints();

	if ( dump_snapshots_ ) {
		++snapshot_counter_;
		if ( snapshot_counter_ % snapshot_interval_ == 0 ) {
			std::string pdb_fn = snapshot_prefix_;
			pdb_fn += "_" + ObjexxFCL::format::I(4,4,snapshot_counter_);
			pdb_fn += ".pdb";
			pose.dump_pdb(pdb_fn);
		}
	}

	pack_full_repack_->apply(pose);
	minimizer_->run( pose, mm_, *scorefxn_, *options_ );

	if ( dump_snapshots_ ) {
		++snapshot_counter_;
		if ( snapshot_counter_ % snapshot_interval_ == 0 ) {
			std::string pdb_fn = snapshot_prefix_;
			pdb_fn += "_" + ObjexxFCL::format::I(4,4,snapshot_counter_);
			pdb_fn += ".pdb";
			pose.dump_pdb(pdb_fn);
		}
	}
}

void BackboneTorsionPerturbation::perturb(core::pose::Pose & pose,
	core::Size level, // level 1 is the base
	core::Real max_delta_torsion,
	core::Size local,
	bool rama_biased,
	bool repack,
	bool minimize) {
	if ( level == 1 ) {
		perturbed_res_ = numeric::random::rg().random_range(1, pose.size());
	}
	for ( core::Size ires=1; ires <= pose.size() ; ++ires ) {
		if ( local != 0 ) {
			if ( std::abs((int) ires - (int) perturbed_res_) > (int) ((local - 1) * level/2) ) continue;
		}

		if ( pose.residue_type(ires).is_protein() ) {
			core::Real phi(0), psi(0);
			if ( level == 1 && rama_biased ) {
				core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

				rama.random_phipsi_from_rama( pose.residue_type(ires).aa(), phi, psi);
				if ( pose.residue(ires).has_property( "D_AA" ) ) {
					phi *= -1.0;
					psi *= -1.0;
				}
			} else {
				phi = (2.* numeric::random::rg().uniform() - 1.) * max_delta_torsion + pose.phi(ires);
				psi = (2.* numeric::random::rg().uniform() - 1.) * max_delta_torsion + pose.psi(ires);
			}
			pose.set_phi(ires, phi);
			pose.set_psi(ires, psi);
		}
	}

	pose.conformation().detect_disulfides();
	if ( repack ) {
		pack_full_repack_->apply(pose);
	}
	if ( minimize ) {
		minimizer_->run( pose, mm_, *scorefxn_, *options_ );
	}
}

void BackboneTorsionPerturbation::apply( core::pose::Pose & pose ) {
	try {
		set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" ) );
	} catch( utility::excn::Exception const & e ) {
		TR.Warning << "Unable to create Talaris2013 legacy scorefunction during initialization of BackboneTorsionPerturbation mover.  Will attempt to use default scorefunction instead.  Error message was:\n" << TR.Red << e.msg() << TR.Reset << std::endl;
		set_scorefunction( core::scoring::get_score_function() );
		TR.Warning << "Recovering and carrying on with default scorefunction." << std::endl;
	}
	//core::pack::task::TaskFactoryOP local_tf = new core::pack::task::TaskFactory();
	//local_tf->push_back(new core::pack::task::operation::RestrictToRepacking());
	// if present, task_factory_ always overrides/regenerates task_
	core::pack::task::PackerTaskOP task;
	if ( task_factory_ != nullptr ) {
		task = task_factory_->create_task_and_apply_taskoperations( pose );
		task->restrict_to_repacking();
	} else {
		task = core::pack::task::TaskFactory::create_packer_task( pose );
		task->initialize_from_command_line().restrict_to_repacking();
	}
	pack_full_repack_ = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >( scorefxn_, task );
	//task->show_all_residue_tasks();

	pose.conformation().detect_disulfides();

	// initialize monte carlo
	protocols::moves::MonteCarloOP mc;
	core::Size ncycles;
	core::Size counters(0);

	perturbed_res_ = 0;
	mc = utility::pointer::make_shared< protocols::moves::MonteCarlo >( pose, *scorefxn_, temperature_ );
	mc->set_autotemp( false, temperature_ );
	ncycles = pose.size() * increase_cycles_;

	minimizer_ = utility::pointer::make_shared< core::optimization::AtomTreeMinimizer >();
	//minimizer_ = new core::optimization::CartesianMinimizer();
	options_ = utility::pointer::make_shared< core::optimization::MinimizerOptions >( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );

	options_->max_iter(20);
	mm_.set_bb  ( true );
	mm_.set_chi ( true );
	mm_.set_jump( false );

	if ( dump_snapshots_ ) {
		std::string pdb_fn = snapshot_prefix_;
		pdb_fn += "_" + ObjexxFCL::format::I(4,4,counters);
		pdb_fn += ".pdb";
		pose.dump_pdb(pdb_fn);
	}

	snapshot_counter_ = 0;

	while ( true ) {
		++counters;

		// lowest level move
		perturb(pose, 60.);
		core::Real score=(*scorefxn_)(pose);

		if ( native_ && native_->size() ) {
			core::Real gdtmm(0.);
			core::sequence::SequenceAlignmentOP native_aln;
			gdtmm = get_gdtmm(*native_, pose, native_aln);
			using namespace ObjexxFCL::format;
			TR << "Trace: " << F(8,3,gdtmm) << " " << F(11,3,score) << std::endl;
		}

		mc->boltzmann(pose, "BackboneTorsionPerturbation");

		if ( dump_snapshots_ ) {
			++snapshot_counter_;
			if ( snapshot_counter_ % snapshot_interval_ == 0 ) {
				std::string pdb_fn = snapshot_prefix_;
				pdb_fn += "_" + ObjexxFCL::format::I(4,4,snapshot_counter_);
				pdb_fn += ".pdb";
				pose.dump_pdb(pdb_fn);
			}
		}
		if ( counters >= ncycles ) break;
	}

	mc->show_scores();
	mc->show_counters();

	if ( recover_low_ ) mc->recover_low(pose);
	if ( dump_snapshots_ ) {
		std::string pdb_fn = snapshot_prefix_;
		pdb_fn += "_" + ObjexxFCL::format::I(4,4,counters);
		pdb_fn += ".pdb";
		pose.dump_pdb(pdb_fn);
	}
}

void BackboneTorsionPerturbation::task_factory( core::pack::task::TaskFactoryCOP tf )
{
	runtime_assert( tf != nullptr );
	task_factory_ = tf;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
BackboneTorsionPerturbation::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	if ( tag->hasOption( "start_res" ) ) {
		TR.Warning << "The 'start_res' option of BackboneTorsionPerturbation doesn't actually do anything." << std::endl;
	}
	if ( tag->hasOption( "stop_res" ) ) {
		TR.Warning << "The 'stop_res' option of BackboneTorsionPerturbation doesn't actually do anything." << std::endl;
	}

	if ( tag->hasOption( "native") ) {
		native_ = utility::pointer::make_shared< core::pose::Pose >();
		core::import_pose::pose_from_file( *native_, tag->getOption< std::string >( "native" ) , core::import_pose::PDB_file);
	}

	if ( tag->hasOption( "increase_cycles" ) ) increase_cycles_ = tag->getOption< core::Real >( "increase_cycles" );
	if ( tag->hasOption( "recover_low" ) ) recover_low_ = tag->getOption< bool >( "recover_low" );
	if ( tag->hasOption( "temp" ) ) temperature_ = tag->getOption< core::Real >( "temp" );
	local_ = tag->getOption< core::Size >( "local" , 0);

	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		set_scorefunction ( (datamap.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone() );
	}

	core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	task_factory( new_task_factory );

	if ( tag->hasOption( "dump_snapshots" ) ) {
		dump_snapshots_ = tag->getOption< core::Size >( "dump_snapshots" );
		snapshot_prefix_ = tag->getOption< std::string >( "snapshot_prefix", "snapshot");
		snapshot_interval_ = tag->getOption< core::Size >( "snapshot_interval" , 100);
	}
}

moves::MoverOP BackboneTorsionPerturbation::clone() const {
	return utility::pointer::make_shared< BackboneTorsionPerturbation >( *this );
}
moves::MoverOP BackboneTorsionPerturbation::fresh_instance() const {
	return utility::pointer::make_shared< BackboneTorsionPerturbation >();
}





std::string BackboneTorsionPerturbation::get_name() const {
	return mover_name();
}

std::string BackboneTorsionPerturbation::mover_name() {
	return "BackboneTorsionPerturbation";
}

void BackboneTorsionPerturbation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "start_res", xsct_positive_integer, "nominally lower bound for a region to sample (NOTE: THIS PARAMETER IS IGNORED.)")
		+ XMLSchemaAttribute( "stop_res", xsct_positive_integer, "nominally upper bound for a region to sample (NOTE: THIS PARAMETER IS IGNORED.)")
		+ XMLSchemaAttribute( "native", xs_string, "file path to native pose, used to calculate gdtmm")
		+ XMLSchemaAttribute( "increase_cycles", xsct_real, "multiply cycle count by this")
		+ XMLSchemaAttribute( "recover_low", xsct_rosetta_bool, "recover the lowest-energy structure seen in the Monte Carlo trajectory at the end?")
		+ XMLSchemaAttribute( "temp", xs_decimal, "Monte Carlo temperature")
		+ XMLSchemaAttribute( "scorefxn", xs_string, "use this scorefunction (from the SCOREFXN section)");

	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	std::string const dump_snapshots_warning(" Note that dump_snapshots must be set to make snapshot_prefix or snapshot_interval active");

	attlist
		+ XMLSchemaAttribute( "dump_snapshots", xsct_non_negative_integer, "dump structures during sampling," + dump_snapshots_warning)
		+ XMLSchemaAttribute::attribute_w_default( "snapshot_prefix", xs_string, "prefix for structures dumped during sampling," + dump_snapshots_warning, "snapshot")
		+ XMLSchemaAttribute::attribute_w_default( "snapshot_interval", xsct_positive_integer, "how frequently to dump structures during sampling," + dump_snapshots_warning, "100");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Not documented.  Appears to be some sort of relax analogue.  Related to BackboneTorsionSampler", attlist );
}

std::string BackboneTorsionPerturbationCreator::keyname() const {
	return BackboneTorsionPerturbation::mover_name();
}

protocols::moves::MoverOP
BackboneTorsionPerturbationCreator::create_mover() const {
	return utility::pointer::make_shared< BackboneTorsionPerturbation >();
}

void BackboneTorsionPerturbationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BackboneTorsionPerturbation::provide_xml_schema( xsd );
}


} // moves
} // protocols
