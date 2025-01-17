// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/SilentFileFullModelInputter.hh
/// @brief  header file for SilentFileFullModelInputter class
/// @author James Thompson


#ifndef INCLUDED_protocols_jd3_full_model_inputters_SilentFileFullModelInputter_hh
#define INCLUDED_protocols_jd3_full_model_inputters_SilentFileFullModelInputter_hh

//unit headers
#include <protocols/jd3/full_model_inputters/FullModelInputter.hh>

// package headers
#include <protocols/jd3/full_model_inputters/SilentFileFullModelInputter.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileOptions.fwd.hh>

//utility headers
#include <utility/file/FileName.fwd.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

#include <string>

#include <utility/tag/XMLSchemaGeneration.fwd.hh> // AUTO IWYU For XMLSchemaDefinition

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

class SilentFileFullModelInputter : public FullModelInputter
{
public:

	SilentFileFullModelInputter();

	~SilentFileFullModelInputter() override;

	bool job_available_on_command_line() const override;

	/// @brief Constructs a list of FullModelInputSource objects reading from the
	/// -s or -l command line flags. This stores the names of the PDBs that
	/// are to be read in, and it initializes the input tags based on the pdb
	/// names, stripping the path and the extension from the file name.
	FullModelInputSources full_model_input_sources_from_command_line() override;

	FullModelInputSources full_model_input_sources_from_tag(
		utility::options::OptionCollection const & opts,
		utility::tag::TagCOP tag
	) override;

	/// @brief this function is responsible for filling the pose reference with
	/// the pose indicated by the job.  The Job object (within its InnerJob)
	/// contains a PoseCOP.  This function needs to either fill the pose
	/// reference from the InnerJob or, on first demand of a pose from that
	/// InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill
	/// the reference.  This implementation uses pose_from_pdb
	//virtual void pose_from_job( core::pose::Pose & pose, JobOP job );
	core::pose::PoseOP
	full_model_from_input_source(
		FullModelInputSource const & input_source,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) override;

	/// @brief this function returns the SilentStruct that belongs to the given job
	//virtual core::io::silent::SilentStruct const& struct_from_job( JobOP job );

	core::io::silent::SilentFileData const &
	silent_file_data() const { return *sfd_; };

	/// @brief returns the name for the element that will be used in a job-definition
	/// file for a structure originating from a silent file: "Silent"
	static std::string keyname();

	/// @brief returns the schema for the PDB element used in a job-definition file
	/// including all options that govern how a silent struct is loaded.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static void list_options_read( utility::options::OptionKeyList & read_options );

private:

	void
	initialize_sfd_from_options_and_tag(
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag
	);

	void
	initialize_sfd_from_files_and_tags(
		utility::vector1< utility::file::FileName > const & silent_files,
		utility::vector1< std::string > const & tags
	);


private:
	core::io::silent::SilentFileOptionsOP sf_opts_;
	core::io::silent::SilentFileDataOP sfd_;

}; // SilentFileFullModelInputter

} // namespace full_model_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_SilentFileFullModelInputter_HH
