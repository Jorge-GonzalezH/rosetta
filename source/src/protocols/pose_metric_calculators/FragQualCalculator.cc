// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/PoseMetricCalculators/FragQualCalculator.cc
/// @brief calculate rmsd fragment quality given a pose
/// Roughly, fragment quality is number of fragments which are close to a pose in rmsd
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/pose_metric_calculators/FragQualCalculator.hh>

// Project Headers

#include <basic/MetricValue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/Frame.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

//// C++ headers
#include <ObjexxFCL/format.hh>

#include <core/pose/ResidueIndexDescription.hh>
#include <core/pose/selection.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/chemical/ChemicalManager.fwd.hh>


static basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.FragQualCalculator" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_metric_calculators {


/// @brief default constructor
FragQualCalculator::FragQualCalculator():
	rmsd_cutoff_goodfrag_( 1.0 ),
	ratio_cutoff_goodfrag_( 0.3 ),
	total_goodfrags_( 0 ),
	coverage_( 0 ),
	frag_( /* NULL */ ),
	begin_( nullptr ),
	end_( nullptr ),
	verbose_( false )
{
	goodfrags_.clear();
}

/// @brief value constructor
FragQualCalculator::FragQualCalculator(
	FragSetOP frag,
	Real const rmsd,
	Real const ratio ):
	rmsd_cutoff_goodfrag_( rmsd ),
	ratio_cutoff_goodfrag_( ratio ),
	total_goodfrags_( 0 ),
	coverage_( 0 ),
	frag_(std::move( frag )),
	begin_( nullptr ),
	end_( nullptr ),
	verbose_( false )
{
	goodfrags_.clear();
}

/// @brief copy constructor
FragQualCalculator::FragQualCalculator( FragQualCalculator const & ) = default;

/// @brief destructor
FragQualCalculator::~FragQualCalculator()= default;


/// @brief set fragments
void
FragQualCalculator::set_fragset( FragSetOP const & frag )
{
	frag_ = frag;
}

/// @brief rmsd cutoff of good fragments
void
FragQualCalculator::rmsd_cutoff( Real const & val )
{
	rmsd_cutoff_goodfrag_ = val;
}

/// @brief
void
FragQualCalculator::ratio_cutoff( Real const & val )
{
	ratio_cutoff_goodfrag_ = val;
}

/// @brief
void
FragQualCalculator::set_region( core::Size const val1, core::Size const val2 )
{
	begin_ = core::pose::make_rid_posenum( val1 );
	end_ = core::pose::make_rid_posenum( val2 );
}

void
FragQualCalculator::begin( core::Size const begin )
{
	begin_ = core::pose::make_rid_posenum( begin );
}

void
FragQualCalculator::end( core::Size const end )
{
	end_ = core::pose::make_rid_posenum( end );
}

/// @brief
void
FragQualCalculator::lookup(
	String const & key,
	MetricValueBase * valptr
) const
{
	using namespace core;
	if ( key == "num_goodfrag" ) {
		basic::check_cast( valptr, &total_goodfrags_, "number of fragments within rmsd cutoff " );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_goodfrags_ );
	} else if ( key == "coverage" ) {
		basic::check_cast( valptr, &coverage_, "ratio of the region where good fragments are included more than XXX% " );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( coverage_ );
	} else {
		TR.Fatal << "FragQualCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


// @brief
std::string
FragQualCalculator::print( String const & key ) const
{
	String result;
	if ( key == "num_goodfrag" ) {
		result = utility::to_string( total_goodfrags_ );
	} else if ( key == "coverage" ) {
		result = utility::to_string( coverage_ );
	} else {
		basic::Error() << "FragQualCalculator cannot compute metric " << key << std::endl;
	}
	return result;
} // apply


/// @brief recomute ncontacts
void
FragQualCalculator::recompute( Pose const & pose )
{
	using ObjexxFCL::format::RJ;
	using ObjexxFCL::format::F;
	using core::scoring::CA_rmsd;
	using core::fragment::ConstFrameIterator;

	// initialization
	core::Size begin( 1 ), end( pose.size() );
	if ( begin_ != nullptr ) begin = begin_->resolve_index( pose );
	if ( end_ != nullptr ) end = end_->resolve_index( pose );
	total_goodfrags_ = 0;
	coverage_ = 0;

	goodfrags_.resize( pose.size() );
	for ( core::Size i=1; i<=pose.size(); i++ ) {
		goodfrags_[ i ] = 0;
	}
	utility::vector1< bool > frag_region( pose.size(), false );
	utility::vector1< bool > is_covered( pose.size(), false );

	Pose input_pose( pose ), test_pose( pose );
	core::util::switch_to_residue_type_set( input_pose, core::chemical::CENTROID_t );
	core::util::switch_to_residue_type_set(  test_pose, core::chemical::CENTROID_t );

	for ( ConstFrameIterator frame = frag_->begin(); frame != frag_->end(); ++frame ) {

		core::Size const start ( frame->start() );
		runtime_assert( start <= pose.size() );

		if ( begin > start || end < start ) continue;
		frag_region[ start ] = true;

		for ( core::Size i=1; i<=frame->nr_frags(); i++ ) {
			// insert fragment
			frame->apply( i, test_pose );
			// calc rmsd
			Real rmsd = CA_rmsd( input_pose, test_pose, start, start + frame->length() - 1 );
			if ( rmsd <= rmsd_cutoff_goodfrag_ ) {
				goodfrags_[ start ] += 1;
			}
			if ( verbose_ ) {
				TR << "FRAG_SCORE " << RJ( 6, frame->length() ) << RJ( 6, frame->start() ) << RJ( 6, i ) << F( 10, 4, rmsd ) << std::endl;
			}
		}

		if ( goodfrags_[ start ] >= frame->nr_frags()*ratio_cutoff_goodfrag_ ) {
			is_covered[ start ] = true;
		}
	} // ConstFrameIterator

	// calc coverage
	core::Size count( 0 );
	for ( core::Size i=1; i<=pose.size(); i++ ) {

		if ( frag_region[ i ] ) {
			count ++;
			if ( is_covered[ i ] ) {
				coverage_++;
			}
		}
		TR.Debug << i << " " << goodfrags_[ i ] << std::endl;
		total_goodfrags_ += Real( goodfrags_[ i ] );
	}

	coverage_ = coverage_ / count;

}

/// @brief parse xml
void
FragQualCalculator::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	rmsd_cutoff_goodfrag_ = tag->getOption<Real>( "rmsd_cutoff", 1.0 );
	ratio_cutoff_goodfrag_ = tag->getOption<Real>( "ratio_cutoff", 0.3 );

	begin_ = core::pose::parse_resnum( tag->getOption<std::string>( "begin", "1" ) );
	if ( tag->hasOption( "end" ) ) {
		end_ = core::pose::parse_resnum( tag->getOption<std::string>( "end" ) );
	} else {
		end_ = utility::pointer::make_shared< core::pose::ResidueIndexDescriptionLastResidue >();
	}

	String const fset_string ( tag->getOption<String>( "frag", "" ) );
	runtime_assert( ! fset_string.empty() );
	if ( data.has( "fragsets", fset_string ) ) {
		frag_ = data.get_ptr<FragSet>( "fragsets", fset_string );
	} else {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "fragsets " + fset_string + " not found in basic::datacache::DataMap.");
	}

	verbose_ = tag->getOption<bool>( "verbose", false );
}


} // filters
} // protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pose_metric_calculators::FragQualCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( rmsd_cutoff_goodfrag_ ) ); // Real
	arc( CEREAL_NVP( ratio_cutoff_goodfrag_ ) ); // Real
	arc( CEREAL_NVP( total_goodfrags_ ) ); // Real
	arc( CEREAL_NVP( coverage_ ) ); // Real
	arc( CEREAL_NVP( goodfrags_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( frag_ ) ); // FragSetOP
	arc( CEREAL_NVP( begin_ ) ); // ResidueIndexDescriptionCOP
	arc( CEREAL_NVP( end_ ) ); // ResidueIndexDescriptionCOP
	arc( CEREAL_NVP( verbose_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pose_metric_calculators::FragQualCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( rmsd_cutoff_goodfrag_ ); // Real
	arc( ratio_cutoff_goodfrag_ ); // Real
	arc( total_goodfrags_ ); // Real
	arc( coverage_ ); // Real
	arc( goodfrags_ ); // utility::vector1<core::Size>
	arc( frag_ ); // FragSetOP
	core::pose::ResidueIndexDescriptionOP begin, end;
	arc( begin ); // core::Size
	begin_ = begin;
	arc( end ); // core::Size
	end_ = end;
	arc( verbose_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pose_metric_calculators::FragQualCalculator );
CEREAL_REGISTER_TYPE( protocols::pose_metric_calculators::FragQualCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_FragQualCalculator )
#endif // SERIALIZATION
