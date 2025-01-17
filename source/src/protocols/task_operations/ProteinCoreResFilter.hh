// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictResiduesToRepackingOperation.hh
/// @author Gabi Pszolla & Sarel Fleishman

#ifndef INCLUDED_protocols_task_operations_ProteinCoreResFilter_hh
#define INCLUDED_protocols_task_operations_ProteinCoreResFilter_hh

// Unit Headers

#include <core/pack/task/operation/ResFilter.hh>

// Project Headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <core/types.hh>

#include <string>

namespace protocols {
namespace task_operations {

class ProteinCore : public core::pack::task::operation::ResFilter {
public:
	ProteinCore();
	bool operator() ( Pose const &, core::Size ) const override;
	core::pack::task::operation::ResFilterOP clone() const override {return (utility::pointer::make_shared< ProteinCore >( *this));}
	void parse_tag( TagCOP ) override;

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real distance_threshold_; // dflt 8.0A; sphere around the residue
	core::Size neighbor_cutoff_; // dflt 10; how many residues in the sequence around the target residue to ignore in computing neighbours
	bool bound_; // dflt false; treat the bound or unbound pose
	core::Size jump_; // dflt 1;
	core::Size neighbor_count_cutoff_; //dflt 6
};
} //namespace protocols
} //namespace task_operations

#endif // INCLUDED_protocols_TaskOperations_ProteinCoreResFilter_HH
