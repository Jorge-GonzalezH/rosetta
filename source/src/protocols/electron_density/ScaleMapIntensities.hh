// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file

#ifndef INCLUDED_protocols_electron_density_ScaleMapIntensities_hh
#define INCLUDED_protocols_electron_density_ScaleMapIntensities_hh

#include <protocols/moves/Mover.hh>



namespace protocols {
namespace electron_density {


/// scale density map intensities to match a pose's
class ScaleMapIntensities : public moves::Mover {
public:
	ScaleMapIntensities();
	~ScaleMapIntensities() override = default;

	void init();

	void apply( core::pose::Pose & ) override;


	moves::MoverOP clone() const override { return utility::pointer::make_shared< ScaleMapIntensities >( *this ); }
	moves::MoverOP fresh_instance() const override { return utility::pointer::make_shared< ScaleMapIntensities >(); }

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real res_low_, res_high_, fade_width_, b_sharpen_;
	core::Size nresbins_;
	bool asymm_only_, ignore_bs_, bin_squared_, mask_, mask_output_, truncate_only_;
	std::string outmap_name_;
};

} // moves
} // protocols

#endif
