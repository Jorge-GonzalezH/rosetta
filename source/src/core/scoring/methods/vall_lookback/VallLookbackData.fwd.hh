// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  
/// @brief
/// @author TJ Brunette


#ifndef INCLUDED_core_scoring_methods_VallLookbackData_fwd_hh
#define INCLUDED_core_scoring_methods_VallLookbackData_fwd_hh

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace methods {

class VallLookbackData;
typedef utility::pointer::shared_ptr< VallLookbackData > VallLookbackDataOP;
typedef utility::pointer::shared_ptr< VallLookbackData const > VallLookbackDataCOP;

}//end methods
}//end scoring
}//end core
#endif
