// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_RepScreenInfo.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_Grid.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/farna/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>



static basic::Tracer TR( "protocols.stepwise.modeler.rna.checker.VDW_CachedRepScreenInfo" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


// @brief Constructor
VDW_CachedRepScreenInfo::VDW_CachedRepScreenInfo() :
	CacheableData(),
	VDW_screen_bin_( VDW_GridCOP( new VDW_Grid() ) )
{
	read_in_VDW_rep_screen_pose_from_command_line();
}

// @brief Constructor
VDW_CachedRepScreenInfo::VDW_CachedRepScreenInfo( core::pose::Pose const & pose ) :
	CacheableData()
{
	VDW_rep_screen_info_list_ = const_vdw_cached_rep_screen_info_from_pose( pose ).VDW_rep_screen_info_list();
	VDW_screen_bin_ = const_vdw_cached_rep_screen_info_from_pose( pose ).VDW_screen_bin();
}

/// @details Copy constructors must copy all data, not just some...
VDW_CachedRepScreenInfo::VDW_CachedRepScreenInfo( VDW_CachedRepScreenInfo const & src ) :
	CacheableData(),
	VDW_rep_screen_info_list_( src.VDW_rep_screen_info_list_ ),
	VDW_screen_bin_( src.VDW_screen_bin_ ) // note, we are not cloning -- this is 'scratch' space that takes a while to initialize.
{
	// we will not copy over VDW_screen_bin, which is a huge 3D grid, and
	// most of the time is zero.
	// only instantiate that when we need it.
}


// @brief default destructor
VDW_CachedRepScreenInfo::~VDW_CachedRepScreenInfo()
{}


basic::datacache::CacheableDataOP
VDW_CachedRepScreenInfo::clone() const
{
	return basic::datacache::CacheableDataOP( new VDW_CachedRepScreenInfo( *this ) );
}


utility::vector1< VDW_RepScreenInfo > &
VDW_CachedRepScreenInfo::VDW_rep_screen_info_list() const
{
	return VDW_rep_screen_info_list_;
}


VDW_GridCOP
VDW_CachedRepScreenInfo::VDW_screen_bin() const
{
	return VDW_screen_bin_;
}


void
VDW_CachedRepScreenInfo::read_in_VDW_rep_screen_pose( VDW_RepScreenInfo & VDW_rep_screen_info ) const
{
	using namespace core::chemical;

	TR.Debug << "importing VDW_rep_screen_pose: " << VDW_rep_screen_info.pose_name << std::endl;
	ResidueTypeSetCOP rsd_set( /*core::chemical::*/ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	if ( VDW_rep_screen_info.pose_name == "" ) utility_exit_with_message( VDW_rep_screen_info.pose_name == "" );
	VDW_rep_screen_info.VDW_pose = core::pose::PoseOP( new core::pose::Pose );

	import_pose::pose_from_pdb( *VDW_rep_screen_info.VDW_pose, *rsd_set, VDW_rep_screen_info.pose_name );
	protocols::farna::make_phosphate_nomenclature_matches_mini( *VDW_rep_screen_info.VDW_pose );
	add_virtual_O2Prime_hydrogen( *VDW_rep_screen_info.VDW_pose );
}


void
VDW_CachedRepScreenInfo::read_in_VDW_rep_screen_pose_from_command_line() const
{
	using namespace basic::options;

	if ( !option[ OptionKeys::stepwise::rna::VDW_rep_screen_info ].user() ) return;
	utility::vector1< std::string > const & All_VDW_rep_screen_pose_info = option[ OptionKeys::stepwise::rna::VDW_rep_screen_info ]();

	for ( Size n = 1; n <= All_VDW_rep_screen_pose_info.size(); n++ ) {

		if ( All_VDW_rep_screen_pose_info[n].find(".pdb", All_VDW_rep_screen_pose_info[n].size()-4) == std::string::npos ) continue; // must have '.pdb'

		VDW_RepScreenInfo VDW_rep_screen_info = VDW_RepScreenInfo();
		VDW_rep_screen_info.pose_name = All_VDW_rep_screen_pose_info[n];
		read_in_VDW_rep_screen_pose( VDW_rep_screen_info );

		VDW_rep_screen_info_list_.push_back( VDW_rep_screen_info );
	}
}



/// @details Pose must already contain a vdw_cached_rep_screen_info object or this method will fail.
VDW_CachedRepScreenInfo const &
const_vdw_cached_rep_screen_info_from_pose( core::pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) );
	return *( utility::pointer::static_pointer_cast< VDW_CachedRepScreenInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) ) );
}


/// @details Either returns a non-const reference to the vdw_cached_rep_screen_info object already stored
/// in the pose, or creates a new vdw_cached_rep_screen_info object, places it in the pose, and returns
/// a non-const reference to it.
VDW_CachedRepScreenInfo &
nonconst_vdw_cached_rep_screen_info_from_pose( core::pose::Pose & pose )
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) ) {
		pose.data().set( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO, VDW_CachedRepScreenInfoOP( new VDW_CachedRepScreenInfo() ) );
	}
	assert( pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) );
	return *( utility::pointer::static_pointer_cast< VDW_CachedRepScreenInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) ) );
}


//////////////////////////////////////////////////////////////////////////////
// basically an alias
VDW_CachedRepScreenInfo const &
make_sure_vdw_cached_rep_screen_info_is_setup( core::pose::Pose & pose )
{
	return nonconst_vdw_cached_rep_screen_info_from_pose( pose );
}


//////////////////////////////////////////////////////////////////////////////
bool
vdw_cached_rep_screen_info_is_setup( core::pose::Pose const & pose )
{
	return pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO );
}


//////////////////////////////////////////////////////////////////////////////
bool
option_vdw_rep_screen_info_user()
{
	return basic::options::option[ basic::options::OptionKeys::stepwise::rna::VDW_rep_screen_info ].user();
}


//////////////////////////////////////////////////////////////////////////////
void
set_vdw_cached_rep_screen_info( core::pose::Pose & pose, VDW_CachedRepScreenInfoOP & vdw_cached_rep_screen_info ){
	pose.data().set( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO, vdw_cached_rep_screen_info );
}


//////////////////////////////////////////////////////////////////////////////
void
set_vdw_cached_rep_screen_info_from_pose( core::pose::Pose & new_pose, core::pose::Pose const & pose ){
	if ( !option_vdw_rep_screen_info_user() ) return;
	VDW_CachedRepScreenInfoOP vdw_cached_rep_screen_info( new VDW_CachedRepScreenInfo( pose ) );
	set_vdw_cached_rep_screen_info( new_pose, vdw_cached_rep_screen_info );
}


///////////////////////////////////////////////////////////////////////////////////////
void
fill_vdw_cached_rep_screen_info_from_command_line( core::pose::Pose & pose ) {
	if ( !option_vdw_rep_screen_info_user() ) return;
	make_sure_vdw_cached_rep_screen_info_is_setup( pose );
}


///////////////////////////////////////////////////////////////////////////////////////
void
fill_vdw_cached_rep_screen_info_from_command_line( utility::vector1< core::pose::Pose * > & pose_pointers ) {
	if ( !option_vdw_rep_screen_info_user() ) return;
	for ( core::Size n = 1; n <= pose_pointers.size(); n++ ) {
		fill_vdw_cached_rep_screen_info_from_command_line( *pose_pointers[n] );
	}
}


} //checker
} //rna
} //modeler
} //stepwise
} //protocols