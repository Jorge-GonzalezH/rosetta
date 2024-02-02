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
#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

using namespace core::scoring;

utility::vector1< std::pair< core::Size, core::Size > >
		identify_secondary_structure_spans( std::string const & ss_string ){
              
     utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
     core::Size strand_start = -1;
     for (core::Size ii = 0; ii < ss_string.size(); ++ii) {
        if (ss_string[ii] == 'E' || ss_string[ii] == 'H') {
          if (int(strand_start) == -1) {
            strand_start = ii;
          } else if (ss_string[ii] != ss_string[strand_start]) {
             ss_boundaries.push_back(std::make_pair(strand_start + 1, ii));
             strand_start = ii;
          }
        } else {
          if (int(strand_start) != -1) {
             ss_boundaries.push_back(std::make_pair(strand_start + 1, ii));
             strand_start = -1;
           }
        }
    }
    if (int(strand_start) != -1) {
      //                                                                                                                                
      ss_boundaries.push_back(std::make_pair(strand_start + 1, ss_string.size()));
    }
    
    for (core::Size ii = 1; ii <= ss_boundaries.size(); ++ii) {
       std::cout << "SS Element " << ii << " from residue "
         << ss_boundaries[ii].first << " to "
         << ss_boundaries[ii].second << std::endl;
    	}
    	return ss_boundaries;
		}
		
int main(int argc, char **argv) {
    devel::init(argc, argv);

    utility::vector1<std::string> filenames = basic::options::option[basic::options::OptionKeys::in::file::s].value();
    if (filenames.size() > 0) {
        std::cout << "You entered: " << filenames[1] << " as the PDB file to be read" << std::endl;
    } else {
        std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }
    std::string secstruct_codes = "EHHEHH";
    std::string ss_string = "EEHHEHHE";

    utility::vector1< char > caracter;
    //se agrego en la sección 3 
    for ( std::string::size_type ii = 0; ii < secstruct_codes.size(); ++ii ) {
    		caracter.push_back(secstruct_codes[ii]);
    }
    
    

    ScoreFunctionOP sfxn = get_score_function();
    core::pose::PoseOP mypose = core::import_pose::pose_from_file(filenames[1]);
    core::Real score = sfxn->score(*mypose); // expect an error
    std::cout << score << std::endl;
    
    //Pymol mover
    protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
    the_observer->pymol().apply( *mypose);
    //a- random number generator
    //core::Size randres = static_cast<core::Size>(numeric::random::random_range(1, mypose->total_residue()));
    
    // Initialize a MonteCarlo object
    protocols::moves::MonteCarlo mc(*mypose, *sfxn, 1.0);
    
    
		core::pose::Pose copy_pose;
    // Begin loop
    // Declarar e inicializar los contadores de verdadero y falso
		int contadorTrue = 0;
		int contadorFalse = 0;
		// Declarar e inicializar el contador
    static int call_count = 0; // Static para que persista entre llamadas de la función

    for (int i = 0; i < 100; ++i) { // Adjust the loop iterations as needed
        // Perturb the phi/psi values for your Pose
        core::Size randres = static_cast<core::Size>(numeric::random::random_range(1, mypose->total_residue()));
        core::Real pert1 = static_cast<core::Real>(numeric::random::random_range(1, 180));
        core::Real pert2 = static_cast<core::Real>(numeric::random::random_range(1, 180));
        core::Real orig_phi = mypose->phi(randres);
        core::Real orig_psi = mypose->psi(randres);
        mypose->set_phi(randres, orig_phi + pert1);
        mypose->set_psi(randres, orig_psi + pert2);

        //you will add packing and minimization calls to your application after you have perturbed the phi and psi values.
        
        core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
				repack_task->restrict_to_repacking();
				core::pack::pack_rotamers( *mypose, *sfxn, repack_task );

        //but before you have called the MonteCarlo::boltzmann function.
        core::kinematics::MoveMap mm;
				mm.set_bb( true );
				mm.set_chi( true );
				
        core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
        
        core::optimization::AtomTreeMinimizer atm;
        //atm.run( *mypose, mm, *sfxn, min_opts );
        
        copy_pose = *mypose;        // crear esta copia hace al loop mucho mas rapido
        atm.run( copy_pose, mm, *sfxn, min_opts );
        *mypose = copy_pose;
        
        // Call MonteCarlo object’s boltzmann method, passing it your Pose
        mc.boltzmann(*mypose);
        // Add code to measure the acceptance / rejection counter
        bool is_accepted_ = mc.mc_accepted();
        
        // Incrementar el contador en cada llamada
        ++call_count;
       

				// Incrementar los contadores según si se acepta o no
				if (is_accepted_) {
			    contadorTrue += 1;
				} else {
			    contadorFalse += 1;
				}

				// Bucle while incorrectamente estructurado, debe ser un if para la impresión
				//if (contadorTrue < 100) {
		    std::cout << "True Counter: " << contadorTrue << std::endl;
 			  std::cout << "False Counter: " << contadorFalse << std::endl;

				//}

				// Calcula la tasa de aceptación
				double acceptanceRate = static_cast<double>(contadorTrue) / (contadorFalse + contadorTrue);

				// Imprime la tasa de aceptación
					std::cout << "Acceptance ratio is: " << acceptanceRate << std::endl;

				//return acceptanceRate;
    }


    // Output final score
    core::Real final_score = sfxn->score(*mypose);
    std::cout << "Final score: " << final_score << std::endl;

    return 0;
}

