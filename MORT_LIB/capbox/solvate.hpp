#ifndef MORTSRC_CAPBOX_SOLVATE_HPP
#define MORTSRC_CAPBOX_SOLVATE_HPP

#include <common.hpp>
#include "solute.hpp"

namespace mort
{
    class molecule_t;

    numvec solvatecap_core( solute_i& solute, const molecule_t& m, const numvec& capcnt, double caprad, double closeness );
 
    numvec solvatebox_core( solute_i& solute, const molecule_t& m, double buffer, double closeness );
 
    double solvateoct_core( solute_i& solute, const molecule_t& m, double buffer, double closeness );
 
    numvec solvateshl_core( solute_i& solute, const molecule_t& m, double shlext, double closeness );
 
    void solvatecap( molecule_t& mol, const molecule_t& svt, const numvec& capcnt, double caprad, double closeness );
 
    void solvatebox( molecule_t& mol, const molecule_t& svt, double buffer, double closeness );
 
    void solvateoct( molecule_t& mol, const molecule_t& svt, double buffer, double closeness );
 
    void solvateshl( molecule_t& mol, const molecule_t& svt, double shlext, double closeness );

} // namespace mort

#endif

