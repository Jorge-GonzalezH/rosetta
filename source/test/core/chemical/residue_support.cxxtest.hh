// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/residue_support.cxxtest.hh
/// @brief unit tests for the residue_support file
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/Atom.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/types.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <set>
#include <cmath>

static basic::Tracer TR("core.chemical.restype_support.cxxtest");

class residue_support_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_shortest_path() {
		using namespace core::chemical;

		ResidueTypeSetCOP rsd_types = ChemicalManager::get_instance()->residue_type_set(FA_STANDARD);

		MutableResidueTypeCOP trp( new MutableResidueType( *rsd_types->name_mapOP("TRP") ) );

		utility::vector1< VD > mainchain = core::chemical::mainchain_path( *trp );
		utility::vector1< VD > mainchain_ref{ trp->atom_vertex("N"), trp->atom_vertex("CA"), trp->atom_vertex("C") };
		TS_ASSERT_EQUALS( mainchain, mainchain_ref );

		utility::vector1< VD > path1 = core::chemical::shortest_path( *trp, trp->atom_vertex("CD1"), trp->atom_vertex("CZ3") );
		utility::vector1< VD > path2_ref{ trp->atom_vertex("CD1"), trp->atom_vertex("CG"), trp->atom_vertex("CD2"), trp->atom_vertex("CE3"), trp->atom_vertex("CZ3") };
		TS_ASSERT_EQUALS( path1, path2_ref );

	}

	void test_renaming() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		ResidueTypeSetCOP rsd_types = cm->residue_type_set(tag);

		MutableResidueTypeCOP rsd_ref;
		rsd_ref = MutableResidueTypeOP( new MutableResidueType( *rsd_types->name_mapOP("LYS") ) );
		MutableResidueTypeOP rsd;
		std::set< std::string > names;

		// Already named - should be no changes.
		rsd = utility::pointer::make_shared< MutableResidueType >( *rsd_ref );
		rename_atoms( *rsd, /*preserve=*/true);
		TS_ASSERT_EQUALS( rsd->natoms(), rsd_ref->natoms() ); // No change in atom number
		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			// Atoms should have matched order in each.
			TS_ASSERT_EQUALS( rsd->atom( rsd->all_atoms()[ii] ).name(), rsd_ref->atom( rsd_ref->all_atoms()[ii] ).name() );
		}

		// Force renaming.
		rsd = utility::pointer::make_shared< MutableResidueType >( *rsd_ref );
		names.clear();
		rename_atoms( *rsd, /*preserve=*/false);
		for ( VD atm: rsd->all_atoms() ) {
			TS_ASSERT_EQUALS( names.count( rsd->atom(atm).name() ), 0 ); // Names should be unique
			names.insert( rsd->atom(atm).name() );
		}
		//non-terminal Lysine: 6 carbons, 1 oxygen, 2 nitrogen, 13 hydrogens.
		TS_ASSERT_EQUALS( names.size(), 22 );
		TS_ASSERT_EQUALS( names.count( " O1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " O2 " ), 0 );
		TS_ASSERT_EQUALS( names.count( " N1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " N2 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " C1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " C6 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " C7 " ), 0 );
		TS_ASSERT_EQUALS( names.count( " H5 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H13" ), 1 );
		TS_ASSERT_EQUALS( names.count( " H14" ), 0 );

		// Partial renaming.
		rsd = utility::pointer::make_shared< MutableResidueType >( *rsd_ref );
		//rsd->rename_atom(rsd->atom_vertex(" NZ "), " N  "); // We don't permit duplicate names anymore
		rsd->rename_atom(rsd->atom_vertex("1HB "), "");
		rsd->rename_atom(rsd->atom_vertex("2HB "), "");
		//rsd->rename_atom(rsd->atom_vertex(" HA "), " H  "); // We don't permit duplicate names anymore
		rename_atoms( *rsd, /*preserve=*/true);
		names.clear();
		for ( VD atm: rsd->all_atoms() ) {
			TS_ASSERT_EQUALS( names.count( rsd->atom(atm).name() ), 0 ); // Names should be unique
			names.insert( rsd->atom(atm).name() );
		}
		TS_ASSERT_EQUALS( names.size(), 22 );
		TS_ASSERT_EQUALS( names.count( " O  " ), 1 );
		TS_ASSERT_EQUALS( names.count( " CA " ), 1 );
		//TS_ASSERT_EQUALS( names.count( " N  " ), 0 ); // We don't permit duplicate names anymore
		//TS_ASSERT_EQUALS( names.count( " NZ " ), 0 ); // We don't permit duplicate names anymore
		//TS_ASSERT_EQUALS( names.count( " N1 " ), 1 ); // We don't permit duplicate names anymore
		//TS_ASSERT_EQUALS( names.count( " N2 " ), 1 ); // We don't permit duplicate names anymore
		TS_ASSERT_EQUALS( names.count( " H1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H2 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H3 " ), 0 );
		//TS_ASSERT_EQUALS( names.count( " H4 " ), 1 );
		//TS_ASSERT_EQUALS( names.count( " H5 " ), 0 );
		//TS_ASSERT_EQUALS( names.count( " H  " ), 0 );
		//TS_ASSERT_EQUALS( names.count( " HA " ), 0 );
		TS_ASSERT_EQUALS( names.count( ""     ), 0 );
		TS_ASSERT_EQUALS( names.count( "1HZ " ), 1 );

	}

	void test_rigid_matrix() {
		using namespace core;
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		// minirosett_database/chemical/atom_type_sets/<tag>
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		// minirosetta_database/chemical/element_sets/<tag>
		ElementSetCOP element_types = cm->element_set("default");
		// minirosetta_database/chemical/mm_atom_type_sets/<tag>
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);

		MutableResidueTypeOP rsd( new MutableResidueType( atom_types, element_types, mm_atom_types, NULL ) );

		VD vd1 = rsd->add_atom( "C1", "aroC", "VIRT", 0 );
		VD vd2 = rsd->add_atom( "C2", "aroC", "VIRT", 0 );
		VD vd3 = rsd->add_atom( "C3", "aroC", "VIRT", 0 );
		VD vd4 = rsd->add_atom( "C4", "aroC", "VIRT", 0 );
		VD vd5 = rsd->add_atom( "C5", "aroC", "VIRT", 0 );
		VD vd6 = rsd->add_atom( "C6", "aroC", "VIRT", 0 );
		VD vd7 = rsd->add_atom( "C7", "aroC", "VIRT", 0 );
		VD vd8 = rsd->add_atom( "C8", "aroC", "VIRT", 0 );
		VD vd9 = rsd->add_atom( "C9", "aroC", "VIRT", 0 );
		VD vd10 = rsd->add_atom( "H10", "Haro", "VIRT", 0 );
		VD vd11 = rsd->add_atom( "H11", "Haro", "VIRT", 0 );
		VD vd12 = rsd->add_atom( "H12", "Haro", "VIRT", 0 );
		/*VD vd13 =*/ rsd->add_atom( "O13", "OOC", "VIRT", 0 );
		// shape :
		//    8-9-10
		//    #
		// 11 7=6
		//  |   |
		//  2-4-5-13
		//  | |   |
		//  1-3   12
		rsd->add_bond( "C1", "C2" );
		rsd->add_bond( "C1", "C3" );
		rsd->add_bond( "C3", "C4" );
		rsd->add_bond( "C4", "C2" );
		rsd->add_bond( "C4", "C5" );
		rsd->add_bond( "C5", "C6" );
		rsd->add_bond( "C6", "C7", DoubleBond );
		rsd->add_bond( "C7", "C8", TripleBond );
		rsd->add_bond( "C8", "C9" );
		rsd->add_bond( "C9", "H10" );
		rsd->add_bond( "C2", "H11" );
		rsd->add_bond( "C5", "O13" );
		rsd->add_bond( "O13", "H12" );

		rsd->atom("C1").ideal_xyz( Vector(0,0,0) );
		rsd->atom("C2").ideal_xyz( Vector(0,1,0) );
		rsd->atom("C3").ideal_xyz( Vector(1,0,0) );
		rsd->atom("C4").ideal_xyz( Vector(1,1,0) );
		rsd->atom("C5").ideal_xyz( Vector(2.1,1,0) );
		rsd->atom("C6").ideal_xyz( Vector(2.1,2.1,0) );
		rsd->atom("C7").ideal_xyz( Vector(1,2.1,0) );
		rsd->atom("C8").ideal_xyz( Vector(1,3.3,0) );
		rsd->atom("C9").ideal_xyz( Vector(2.5,3.3,0) );
		rsd->atom("H10").ideal_xyz( Vector(3.0,3.3,0) );
		rsd->atom("H11").ideal_xyz( Vector(0,1.5,0) );
		rsd->atom("H12").ideal_xyz( Vector(3.3,0,0) );
		rsd->atom("O13").ideal_xyz( Vector(3.3,1,0) );

		rsd->bond("C1","C2").ringness(BondInRing);
		rsd->bond("C1","C3").ringness(BondInRing);
		rsd->bond("C3","C4").ringness(BondInRing);
		rsd->bond("C2","C4").ringness(BondInRing);

		core::chemical::VDDistanceMatrix distances(*rsd, 1e9);

		calculate_rigid_matrix( *rsd, distances );

		core::Real const delta( 0.001 );
		// Distances should be symetrical.
		for ( VD atm1: rsd->all_atoms() ) {
			for ( VD atm2: rsd->all_atoms() ) {
				TS_ASSERT_DELTA( distances(atm1,atm2), distances(atm2,atm1), delta );
			}
		}

		TS_ASSERT_DELTA( distances(vd1,vd2), 1.0, delta );
		TS_ASSERT_DELTA( distances(vd1,vd4), sqrt( 1.0*1.0 + 1.0*1.0), delta );
		TS_ASSERT_DELTA( distances(vd1,vd5), sqrt(2) + 1.1, delta );
		TS_ASSERT_DELTA( distances(vd5,vd3), 1 + 1.1, delta );
		TS_ASSERT_DELTA( distances(vd1,vd6), sqrt(2) + 1.1 + 1.1, delta );
		TS_ASSERT_DELTA( distances(vd6,vd4), 1.1 + 1.1, delta );
		//Bonding issues
		TS_ASSERT_DELTA( distances(vd7,vd8),1.2, delta );
		TS_ASSERT_DELTA( distances(vd7,vd6),1.1, delta );
		TS_ASSERT_DELTA( distances(vd6,vd8), sqrt( 1.1*1.1 + 1.2*1.2 ), delta );
		TS_ASSERT_DELTA( distances(vd4,vd8), 1.1 + 1.1 + sqrt( 1.1*1.1 + 1.2*1.2 ), delta );
		//C9 is a Non-rotameric stub, so should be non-rotameric
		TS_ASSERT_DELTA( distances(vd9,vd6),sqrt( 0.4*0.4 + 1.2*1.2), delta );
		//Across everything
		TS_ASSERT_DELTA( distances(vd1,vd9), sqrt(2) + 1.1 + 1.1 + sqrt( 0.4*0.4 + 1.2*1.2 ), delta );
		//Hydrogens add to rigid unit.
		TS_ASSERT_DELTA( distances(vd7,vd10),sqrt(2*2 + 1.2*1.2), delta );
		TS_ASSERT_DELTA( distances(vd3,vd11),sqrt(1*1+1.5*1.5), delta );
		//But polar hydrogens don't
		TS_ASSERT_DELTA( distances(vd5,vd12),1.2+1, delta );
	}


	void test_possible_nbr() {
		using namespace core;
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		// minirosett_database/chemical/atom_type_sets/<tag>
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		// minirosetta_database/chemical/element_sets/<tag>
		ElementSetCOP element_types = cm->element_set("default");
		// minirosetta_database/chemical/mm_atom_type_sets/<tag>
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);

		MutableResidueTypeOP rsd( new  MutableResidueType( atom_types, element_types, mm_atom_types, NULL ) );

		rsd->add_atom( "C1", "aroC", "VIRT", 0 );
		rsd->add_atom( "C2", "aroC", "VIRT", 0 );
		rsd->add_atom( "C3", "aroC", "VIRT", 0 );
		rsd->add_atom( "C4", "aroC", "VIRT", 0 );
		rsd->add_atom( "C5", "aroC", "VIRT", 0 );
		rsd->add_atom( "C6", "aroC", "VIRT", 0 );
		// "E" shape, and all rigid
		rsd->add_bond( "C1", "C2", DoubleBond );
		rsd->add_bond( "C2", "C3", DoubleBond );
		rsd->add_bond( "C3", "C4", DoubleBond );
		rsd->add_bond( "C4", "C5", DoubleBond );
		rsd->add_bond( "C3", "C6", DoubleBond );

		rsd->atom("C1").ideal_xyz( Vector(0,0,0) );
		rsd->atom("C2").ideal_xyz( Vector(2,2,0) );
		rsd->atom("C3").ideal_xyz( Vector(3,1,0) );
		rsd->atom("C4").ideal_xyz( Vector(4,0,0) );
		rsd->atom("C5").ideal_xyz( Vector(2,-2,0) );
		rsd->atom("C6").ideal_xyz( Vector(2,0,0) );

		core::Real maxdist(0);
		VD nbr;

		//  // C6 is omitted because it's singly bonded.
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( nbr, rsd->vd_from_name("C3") );
		//  TS_ASSERT_DELTA( maxdist, 3 * sqrt( 2 ), 1e-6 ); //To C1 & C5
		//
		//  // Doubly bonded H is still omitted.
		//  rsd->add_bond( "C6", "C2" );
		//  rsd->atom("C6").element_type( element_types->element("H") );
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C3");
		//  TS_ASSERT_DELTA( maxdist, 3 * sqrt(2) , 1e-6 ); //To C1 & C5
		//
		//  // But a ring through an H will count
		//  rsd->add_bond( "C6", "C2" );
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C6");
		//  TS_ASSERT_DELTA( maxdist, 2.0 , 1e-6 );
		//
		//
		//  rsd->atom("C4").ideal_xyz( Vector(3,2,0) );
		//  rsd->atom("C5").ideal_xyz( Vector(1,3,0) );
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C2" );
		//  TS_ASSERT_DELTA( maxdist, sqrt( 2.0*2.0+2.0*2.0 ), 1e-6 ); //To C1
		//
		//  ////Hydrogens aren't included in furthest distance calculation,
		//  //// but they still count as bonded groups.
		//  //rsd->atom("C1").element_type( element_types->element("H") );
		//  //nbr = ResidueType::null_vertex;
		//  //maxdist = find_nbr_dist(*rsd, nbr);
		//  //TS_ASSERT_EQUALS( nbr, rsd->vd_from_name("C2") );
		//  //TS_ASSERT_DELTA( maxdist, sqrt( 2.0 ), 1e-6 ); //To C3 & C5 - As hydrogen, C6 = sqrt(8) doesn't count.

		// C6 is omitted because it's singly bonded.
		nbr = MutableResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C3" );
		TS_ASSERT_DELTA( maxdist, sqrt( 1.0*1.0 + 3.0*3.0 ), 1e-6 ); //To C1 & C5

		rsd->add_bond( "C6", "C2" );
		nbr = MutableResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C6" );
		TS_ASSERT_DELTA( maxdist, 2.0 , 1e-6 );

		// Doubly bonded H is still omitted.
		rsd->atom("C6").element_type( element_types->element("H") );
		nbr = MutableResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C3" );
		TS_ASSERT_DELTA( maxdist, sqrt( 1.0*1.0 + 3.0*3.0 ) , 1e-6 ); //To C1 & C5

		rsd->atom("C4").ideal_xyz( Vector(3,2,0) );
		rsd->atom("C5").ideal_xyz( Vector(1,3,0) );
		nbr = MutableResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C2" );
		TS_ASSERT_DELTA( maxdist, sqrt( 2.0*2.0+2.0*2.0 ), 1e-6 ); //To C1
	}

	void test_recharging() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		ResidueTypeSetCOP rsd_types = cm->residue_type_set(tag);

		MutableResidueTypeOP rsd( new MutableResidueType( rsd_types->name_map("TYR") ) );

		TR << "Testing standard charging." << std::endl;

		for ( VD atm: rsd->all_atoms() ) {
			rsd->atom(atm).formal_charge( 0.0 );
			rsd->atom(atm).charge( -1.234 ); // Will be reset
		}

		rosetta_recharge_fullatom( *rsd );

		// 21 atoms in TYR -- Nbb + Cabb + CObb + OCbb + CH2 + aroC * 6 + OH + HNbb + Hapo * 3 + Haro * 4 + Hpol
		//core::Real Y_naive_charge = -0.47 + 0.07 + 0.51 + -0.51 + -0.18 + -0.115*6 + -0.66 + 0.31 + 0.115*4 + 0.095*3 + 0.43; // -0.4450
		core::Real Y_naive_charge = -0.47 + 0.07 + 0.51 + -0.51 + -0.18 + -0.115*4 + -0.66 + 0.31 + 0.115*4 + 0.095*3 + 0.43; // -0.2150; two aroC -> CH0 in beta_nov15
		// skip two CH0 with zero-charge: 21 -> 19
		TS_ASSERT_DELTA( rsd->atom(" N  ").charge(), -0.470 - Y_naive_charge/19, 1e-4 ); //"Nbb "
		TS_ASSERT_DELTA( rsd->atom(" CA ").charge(),  0.070 - Y_naive_charge/19, 1e-4 ); //"CAbb"
		TS_ASSERT_DELTA( rsd->atom(" OH ").charge(), -0.660 - Y_naive_charge/19, 1e-4 ); //"OH  "

		core::Real net_charge(0);
		for ( VD atm: rsd->all_atoms() ) {
			net_charge += rsd->atom(atm).charge();
		}
		TS_ASSERT_DELTA( net_charge, 0, 1e-4 );

		TR << "Testing charging with formal charges" << std::endl;

		rsd->atom(" OH ").formal_charge( -1 );
		rsd->atom(" O  ").formal_charge( -1 );

		rosetta_recharge_fullatom( *rsd );

		TS_ASSERT_DELTA( rsd->atom(" N  ").charge(), -0.470 - Y_naive_charge/19 + -2.0/19, 1e-4 ); //"Nbb "
		TS_ASSERT_DELTA( rsd->atom(" CA ").charge(),  0.070 - Y_naive_charge/19 + -2.0/19, 1e-4 ); //"CAbb"
		TS_ASSERT_DELTA( rsd->atom(" OH ").charge(), -0.660 - Y_naive_charge/19 + -2.0/19, 1e-4 ); //"OH  "

		net_charge = 0;
		for ( VD atm: rsd->all_atoms() ) {
			net_charge += rsd->atom(atm).charge();
		}
		TS_ASSERT_DELTA( net_charge, -2.0, 1e-4 );

		TR << "Testing that virtual atoms aren't recharged and aren't counted during recharge" << std::endl;

		rsd = utility::pointer::make_shared< MutableResidueType >( rsd_types->name_map("PRO") );

		for ( VD atm: rsd->all_atoms() ) {
			rsd->atom(atm).formal_charge( 0.0 );
			rsd->atom(atm).charge( -1.234 ); // Will be reset
		}

		rsd->atom(" CB ").formal_charge( -1 );
		rsd->atom(" CG ").formal_charge( -1 );

		rosetta_recharge_fullatom( *rsd );

		// 14 non-virtual in PRO:  Npro + CAbb + CObb + OCbb + CH2 * 3 + Hapo * 7
		core::Real P_naive_charge = -0.37 + 0.07 + 0.51 + -0.51 + -0.18 * 3 + 0.095 * 7;
		TS_ASSERT_DELTA( rsd->atom(" N  ").charge(), -0.370 - P_naive_charge/14 + -2.0/14, 1e-4 ); //"Nbb "
		TS_ASSERT_DELTA( rsd->atom(" CA ").charge(),  0.070 - P_naive_charge/14 + -2.0/14, 1e-4 ); //"CAbb"
		TS_ASSERT_EQUALS( rsd->atom(" NV ").charge(), 0 ); // VIRT
	}

	void test_make_centroid() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP fa_rts = cm->residue_type_set( FA_STANDARD );
		ResidueTypeSetCOP cen_rts = cm->residue_type_set( CENTROID );

		utility::io::izstream paramslist("core/chemical/params/cen_types_list.txt");
		std::string cenfile, fafile;
		paramslist >> cenfile >> fafile;
		while ( paramslist ) {
			TR << "Comparing converted " << fafile << " with " << cenfile << std::endl;

			core::chemical::MutableResidueTypeOP fa_rsd = read_topology_file("core/chemical/params/"+fafile, fa_rts );
			core::chemical::MutableResidueTypeOP cen_rsd = read_topology_file("core/chemical/params/"+cenfile, cen_rts );

			core::chemical::MutableResidueTypeOP converted = make_centroid( *fa_rsd );

			TS_ASSERT_EQUALS( converted->atom_type_set().name(), cen_rsd->atom_type_set().name() );
			TS_ASSERT_EQUALS( converted->natoms(), cen_rsd->natoms() );

			for ( VD atm: cen_rsd->all_atoms() ) {
				std::string const & name( cen_rsd->atom_name( atm ) );
				TSM_ASSERT_EQUALS( name, converted->atom_type( converted->atom_vertex(name) ).name(), cen_rsd->atom_type( cen_rsd->atom_vertex(name) ).name() );
			}
			paramslist >> cenfile >> fafile;
		}

	}

	// Test to make sure that none of the ResidueTypes in the standard FA_STANDARD set crash when being converted.
	void test_make_centroid_no_crash() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP fa_rts = cm->residue_type_set( FA_STANDARD );
		ResidueTypeSetCOP cen_rts = cm->residue_type_set( CENTROID );

		for ( core::chemical::ResidueTypeCOP fa_rsd: fa_rts->base_residue_types() ) {
			TS_ASSERT_THROWS_NOTHING( make_centroid( *fa_rsd ) );
		}
		for ( core::chemical::ResidueTypeCOP fa_rsd: fa_rts->unpatchable_residue_types() ) {
			TS_ASSERT_THROWS_NOTHING( make_centroid( *fa_rsd ) );
		}

	}

};
