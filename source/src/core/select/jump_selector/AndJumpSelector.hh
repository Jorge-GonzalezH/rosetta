// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/AndJumpSelector.hh
/// @brief  The AndJumpSelector combines logic from multiple JumpSelectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_AndJumpSelector_HH
#define INCLUDED_core_select_jump_selector_AndJumpSelector_HH

// Unit headers
#include <core/select/jump_selector/AndJumpSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <list>

namespace core {
namespace select {
namespace jump_selector {

/// @brief The AndJumpSelector combines the output of multiple JumpSelectors using AND
/// logic, i.e., only jumps selected by ALL contained JumpSelectors will be selected.
/// JumpSelecters can be pulled in from a DataMap, from subtags (for JumpSelectors
/// known to the JumpSelectorFactory) or programmatically through %add_jump_selector.
class AndJumpSelector : public JumpSelector {
public:
	// derived from base class
	AndJumpSelector();

	/// @brief Copy constructor
	///
	AndJumpSelector( AndJumpSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	JumpSelectorOP clone() const override;

	AndJumpSelector( JumpSelectorCOP selector1);
	AndJumpSelector( JumpSelectorCOP selector1, JumpSelectorCOP selector2 );

	~AndJumpSelector() override;

	JumpSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//unit-specific
	/**
	* @brief adds a JumpSelector
	*/
	void add_jump_selector(JumpSelectorCOP selector);

	Size num_selectors() const;

	/**
	* @brief applies newSubset to existingSubset and thereby modifies the latter
	*/
	void apply_and_to_subset(JumpSubset const & newSubset, JumpSubset & existingSubset) const;

private: // data members

	std::list< JumpSelectorCOP > selectors_;

};


} //namespace jump_selector
} //namespace select
} //namespace core

#endif
