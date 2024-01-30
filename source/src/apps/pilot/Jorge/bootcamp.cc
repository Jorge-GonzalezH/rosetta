// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include <iostream>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <numeric/random/random.hh>

using namespace core::scoring;
 
int main( int argc, char ** argv ) {
devel::init( argc, argv );

utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
if ( filenames.size() > 0 ) {
std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
} else {
	std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
	return 1;
}

ScoreFunctionOP sfxn = get_score_function();
core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
core::Real score = sfxn->score( *mypose ); // expect an error
std::cout << score << std::endl;



core::Size randres = static_cast< core::Size > (numeric::random::random_range(1, mypose->total_residue())); 

//residue in the pose   
core::Real pert1 = static_cast< core::Size > (numeric::random::random_range(1,180));        //… code here to get a random number
core::Real pert2 = static_cast< core::Size > (numeric::random::random_range(1,180));        //… code here to get another random number
core::Real orig_phi = mypose->phi( randres );
core::Real orig_psi = mypose->psi( randres );
mypose->set_phi( randres, orig_phi + pert1 );
mypose->set_psi( randres, orig_psi + pert2 );

// Call MonteCarlo object’s boltzmann method, passing it your Pose


//Scoring a Pose
//core::Real top = sfxn->score(pose)(*sfxn_)( pose ); 
//core::Real top = sfxn->(*sfxn)( pose ); 

//Want to see the energies? 
//pose.energies().show( std::cout );
//static_cast< core::Size > ( uniform_random_number * N + 1 )

}


