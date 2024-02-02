// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/PerResidueBfactorMetric.cc
/// @brief new metric based on pLDDT of alphafold Bfactor
/// @author Jorge-GonzalezH (jorge.gonzalezh@usm.cl)

// Unit headers
#include <protocols/bootcamp/PerResidueBfactorMetric.hh>
#include <protocols/bootcamp/PerResidueBfactorMetricCreator.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.bootcamp.PerResidueBfactorMetric" );


namespace protocols {
namespace bootcamp {

using namespace core::select;
using namespace core::select::residue_selector;

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
PerResidueBfactorMetric::PerResidueBfactorMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueBfactorMetric::~PerResidueBfactorMetric(){}

core::simple_metrics::SimpleMetricOP
PerResidueBfactorMetric::clone() const {
	return utility::pointer::make_shared< PerResidueBfactorMetric >( *this );
}

std::string
PerResidueBfactorMetric::name() const {
	return name_static();
}

std::string
PerResidueBfactorMetric::name_static() {
	return "PerResidueBfactorMetric";

}
std::string
PerResidueBfactorMetric::metric() const {

	return "SHORT_NAME_FOR_SCOREFILE_HEADER_DEFAULT";
}

void
PerResidueBfactorMetric::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );


	if (tag->hasOption("bogus_option")){
		return;
	}
}

void
PerResidueBfactorMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("bogus_option", xsct_rosetta_bool, "test bogus option", "false");

	//attributes_for_parse_residue_selector( attlist, "residue_selector",
	//	"Selector specifying residues." );

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		"new metric based on pLDDT of alphafold Bfactor", attlist);
}

std::map< core::Size, core::Real >
PerResidueBfactorMetric::calculate(const core::pose::Pose & pose) const {
	utility::vector1< core::Size > selection1 = selection_positions(get_selector()->apply(pose));

}

/// @brief This simple metric is unpublished.  It returns Jorge-GonzalezH as its author.
void
PerResidueBfactorMetric::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"PerResidueBfactorMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Jorge-GonzalezH",
		"TODO: institution",
		"jorge.gonzalezh@usm.cl",
		"Wrote the PerResidueBfactorMetric."
		)
	);
}

void
PerResidueBfactorMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueBfactorMetric::provide_xml_schema( xsd );
}

std::string
PerResidueBfactorMetricCreator::keyname() const {
	return PerResidueBfactorMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueBfactorMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< PerResidueBfactorMetric >();
}

} //bootcamp
} //protocols


#ifdef    SERIALIZATION



template< class Archive >
void
protocols::bootcamp::PerResidueBfactorMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	//arc( CEREAL_NVP( output_as_pdb_nums_ ) );

}

template< class Archive >
void
protocols::bootcamp::PerResidueBfactorMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	//arc( output_as_pdb_nums_ );


}

SAVE_AND_LOAD_SERIALIZABLE( protocols::bootcamp::PerResidueBfactorMetric );
CEREAL_REGISTER_TYPE( protocols::bootcamp::PerResidueBfactorMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_bootcamp_PerResidueBfactorMetric )
#endif // SERIALIZATION




