// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/ConstraintGeneratorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary ConstraintGenerators
///         from a string --> ConstraintGeneratorCreator map
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Package headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <basic/citation_manager/CitationManager.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

static basic::Tracer TR( "protocols.constraint_generator.ConstraintGeneratorFactory" );

namespace protocols {
namespace constraint_generator {

ConstraintGeneratorFactory::ConstraintGeneratorFactory():
	utility::SingletonBase< ConstraintGeneratorFactory >(),
	creator_map_()
{}

void
ConstraintGeneratorFactory::factory_register( ConstraintGeneratorCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string const err_msg = "Factory Name Conflict: Two or more ConstraintGeneratorCreators registered with the name " + creator->keyname();
		utility_exit_with_message(  err_msg );
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool ConstraintGeneratorFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

ConstraintGeneratorOP
ConstraintGeneratorFactory::new_constraint_generator(
	std::string const & constraint_generator_name,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	if ( ! has_type( constraint_generator_name ) ) {
		std::string err_msg =  "No ConstraintGeneratorCreator with the name '" + constraint_generator_name + "' has been registered with the ConstraintGeneratorFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	auto iter = creator_map_.find( constraint_generator_name );
	ConstraintGeneratorOP new_constraint_generator = iter->second->create_constraint_generator();
	new_constraint_generator->parse_my_tag( tag, datamap );
	return new_constraint_generator;
}


/// @brief Get the XML schema for a given constraint generator.
/// @details Throws an error if the constraint generator is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ConstraintGeneratorFactory::provide_xml_schema(
	std::string const &cst_generator_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	auto iter( creator_map_.find( cst_generator_name ) );
	if ( iter != creator_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: ConstraintGeneratorCreatorOP for " + cst_generator_name + " is NULL!" );
		}
		iter->second->provide_xml_schema( xsd );
	} else {
		TR << "Available constraint generators: ";
		for ( auto const & it : creator_map_ ) {
			TR << it.first << ", ";
		}
		TR << std::endl;
		utility_exit_with_message( cst_generator_name + " is not known to the ConstraintGeneratorFactory. Was it registered via a ConstraintGeneratorRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
	}
}


void
ConstraintGeneratorFactory::define_constraint_generator_xml_schema_group( utility::tag::XMLSchemaDefinition & xsd ) const{
	try{
		utility::tag::define_xml_schema_group(
			creator_map_,
			constraint_generator_xml_schema_group_name(),
			& complex_type_name_for_constraint_generator,
			xsd );
	} catch( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for Constraints from ConstraintFactory; offending class"
			" must call protocols::constraint_generator::complex_type_name_for_constraint when defining"
			" its XML Schema\n" + e.msg() );
	}
}

std::string
ConstraintGeneratorFactory::constraint_generator_xml_schema_group_name(){
	return "constraint_generator";
}

std::string
ConstraintGeneratorFactory::complex_type_name_for_constraint_generator( std::string const & constraint_name ){
	return "constraint_generator_" + constraint_name + "_type";
}

void
ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & constraint_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;

	std::string const citation_string( get_instance()->get_citation_humanreadable( constraint_type ) );

	ct_gen.complex_type_naming_func( & complex_type_name_for_constraint_generator )
		.element_name( constraint_type )
		.description( description + "\n\n" + citation_string )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

/// @brief Get a human-readable listing of the citations for a given filter, by filter name.
/// @details Returns an empty string if there are no citations.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
ConstraintGeneratorFactory::get_citation_humanreadable(
	std::string const & constraint_generator_name
) const {
	using namespace basic::citation_manager;
	CitationCollectionList citations;
	auto const iter = creator_map_.find( constraint_generator_name );
	runtime_assert_string_msg( iter != creator_map_.end(), "Error in ConstraintGeneratorFactory::get_citation_humanreadable(): Did not recognize " + constraint_generator_name + "!" );
	ConstraintGeneratorOP new_constraint_generator = iter->second->create_constraint_generator();
	runtime_assert_string_msg( new_constraint_generator != nullptr, "Error in ConstraintGeneratorFactory::get_citation_humanreadable(): Could not instantiate " + constraint_generator_name + "!" );
	new_constraint_generator->provide_citation_info(citations);
	if ( citations.empty() ) return "";
	std::ostringstream ss;
	ss << "References and author information for the " << constraint_generator_name << " constraint generator:" << std::endl;
	ss << std::endl;
	basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info_from_list_to_stream( citations, ss );
	return ss.str();
}

} //namespace constraint_generator
} //namespace protocols
