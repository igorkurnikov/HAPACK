
#include <common.hpp>
#include <object.hpp>

namespace mort
{

    double eval_each_bond(const morf_t& bond)
    {
        assert( bond.has_d(FORCE) && bond.has_d(EQUIL) );
        double value = mort::dist( atom_1st(bond), atom_2nd(bond) );
        double force = bond.get_d(FORCE);
        double equil = bond.get_d(EQUIL);
        return force * (value-equil) * (value-equil);
    }     

    double eval_each_angl( const morf_t& angl )
    {
        assert( angl.has_d(FORCE) && angl.has_d(EQUIL) );
        double value = mort::angl( atom_1st(angl), atom_2nd(angl), atom_3rd(angl) ) * M_PI/180.0;
        double force = angl.get_d(FORCE);
        double equil = angl.get_d(EQUIL);
        return force * (value-equil) * (value-equil);
    }

    double eval_each_tors(const morf_t& tors)
    {
        assert( tors.has_d(FORCE) && tors.has_d(EQUIL) );
        double value = mort::tors( atom_1st(tors), atom_2nd(tors), atom_3rd(tors), atom_4th(tors) ) * M_PI/180.0;
        double force = tors.get_d(FORCE);
        double equil = tors.get_d(EQUIL);
        return force + force * cos( tors.get_i(PERIOD)*value-equil);
    }
        
    double eval_each_oops(const morf_t& oops)
    {
        assert( oops.has_d(FORCE) && oops.has_d(EQUIL) );
        double value = tors( atom_1st(oops), atom_2nd(oops), atom_3rd(oops), atom_4th(oops) ) * M_PI/180.0;
        double force = oops.get_d(FORCE);
        double equil = oops.get_d(EQUIL);
        return force + force * cos(oops.get_i(PERIOD)*value-equil);
    }

    double eval_bond( const molecule_t& mol )
    {
        return mort::sum( mol.bond_begin(), mol.bond_end(), 0.0, std::ptr_fun( &eval_each_bond ) );
    }
    
    double eval_angl( const molecule_t& mol )
    {
        return mort::sum( mol.angl_begin(), mol.angl_end(), 0.0, std::ptr_fun( &eval_each_angl ) );
    }

    double eval_tors( const molecule_t& mol )
    {
        double etor = 0.0;

        diheiter_t tors = mol.dihe_begin();
        for( ; tors != mol.dihe_end(); ++tors )
        {
            etor += eval_each_tors(*tors);
        }

        return etor;
    }

    double eval_oops( const molecule_t& mol )
    {
        return mort::sum( mol.impr_begin(), mol.impr_end(), 0.0, std::ptr_fun( &eval_each_oops ) );
    }

}


