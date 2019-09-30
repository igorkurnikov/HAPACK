#ifndef MORT_COMMON_HPP
#define MORT_COMMON_HPP

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

#ifndef __GNUC__
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#endif

#include "common/numvec.hpp"
#include "common/pertab.hpp"
#include "common/fortran.hpp"
#include "common/stralgo.hpp"
#include "common/crdalgo.hpp"
#include "common/constant.hpp"
#include "common/funstack.hpp"
#include "common/hashcode.hpp"
#include "common/geometry.hpp"

namespace mort
{
    using std::string;
    
    using std::vector;

    using std::istream;

    using std::ostream;

} // namespace mort


#endif


/// \mainpage the programmer's reference of MORT
/// MORT is grouped by the following modules:
///
/// \subpage common
///
/// \subpage object
///
/// \subpage format
///
/// \subpage capbox
///
/// \subpage atmask
///



