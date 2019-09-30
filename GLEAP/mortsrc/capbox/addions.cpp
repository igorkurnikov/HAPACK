#include <stdexcept>
#include <object.hpp>
#include "shaper.hpp"
#include "octree.hpp"
#include "addions.hpp"

#ifndef M_PI
#include <math.h>
#define M_PI 3.1415926535897932385
#endif

namespace mort
{    
    using namespace capbox;

    bool find_bump_solvent( ionee_i& ionee, numvec& pos, double rion, int& rid )
    {
        std::cout << "Info: looking for bump solvent near: " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        vector<double>& pcord = ionee.getcord();
        vector<double>& pvdwr = ionee.getvdwr();
	vector<string>& ptype = ionee.gettype();

        int minid = -1;
        double mind2 = 999.0;

        int bgn = ionee.solute_natom();
        int end = ionee.natom();
        int solvent_size = (end - bgn)/(ionee.nresd()-ionee.solute_nresd());
        for( int i=bgn; i < end; ++i )
        {
            if( ptype[i][0]=='H' || ptype[i][0]=='h' )
            {
                continue;
            }

            double dx = pcord[3*i  ] - pos[0];
            double dy = pcord[3*i+1] - pos[1];
            double dz = pcord[3*i+2] - pos[2];
            double d2 = dx*dx + dy*dy + dz*dz;
            if( minid==-1 || mind2 > d2 )
            {
                minid = i;
                mind2 = d2;
            }
        }

        if( minid==-1 || sqrt(mind2) > (rion+pvdwr[minid]) )
        {
            // no bump
            std::cout << "Info: no bump, min dist: " << sqrt(mind2) << std::endl;
            return false;
        }
        else
        {
            std::cout << "Info: bump with " << pcord[3*minid] << " " << pcord[3*minid+1] << " " << pcord[3*minid+2] << std::endl;
        }


        rid = ionee.solute_nresd() + (minid - bgn)/solvent_size;
        pos[0] = pcord[3*minid];
        pos[1] = pcord[3*minid+1];
        pos[2] = pcord[3*minid+2];
        return true;
    } 


    void minmaxpos( const vector<double>& pcord, numvec& pmin, numvec& pmax, int natom )
    {
        pmin[0] = pcord[0]; pmin[1] = pcord[1]; pmin[2] = pcord[2];
        pmax[0] = pcord[0]; pmax[1] = pcord[1]; pmax[2] = pcord[2];

        for( int i=1; i < natom; ++i )
        {
            pmin[0] = std::min( pmin[0], pcord[3*i] );
            pmin[1] = std::min( pmin[1], pcord[3*i+1] );
            pmin[2] = std::min( pmin[2], pcord[3*i+2] );
            pmax[0] = std::max( pmax[0], pcord[3*i] );
            pmax[1] = std::max( pmax[1], pcord[3*i+1] );
            pmax[2] = std::max( pmax[2], pcord[3*i+2] );
        }
    }


    void addions( molecule_t& m, const molecule_t& ion, int num_ion, double shlext, double resolution )
    {
        mort_ionee ionee(m);
        addions_core( ionee, ion, num_ion, shlext, resolution );
	m.cleanup();
    }

    void addions_core( ionee_i& ionee, const molecule_t& ion, int num_ion, double shlext, double resolution )
    {
        int natom = ionee.natom();
        int nresd = ionee.nresd();
        int solute_nresd = ionee.solute_nresd();
        int solute_natom = ionee.solute_natom();

        bool nosolvent = (natom==solute_natom);

        if( solute_natom != ionee.natom() )
        {
            std::cout << "Info: solute natom, total natom: " << solute_natom << " " << natom << std::endl;
            std::cout << "      solute nresd, total nresd: " << solute_nresd << " " << nresd << std::endl;
        }

        vector<double>& pcord = ionee.getcord();
        vector<double>& pchrg = ionee.getchrg();
        vector<double>& pvdwr = ionee.getvdwr();


        numvec pmin(3), pmax(3);
        minmaxpos( pcord, pmin, pmax, solute_natom );

        double max_radius = *std::max_element( pvdwr.begin(), pvdwr.end() );
        double ion_radius = get_vdwr( ion.atom_begin()->get_s(TYPE) );
        if( ion.natom() > 1 )
        {
            std::cout << "Info: polyatomic ion, radius set to " << ion_radius*2.5 << " (2.5 times first atom)" << std::endl;
            ion_radius *= 2.5;
        }

	if( std::abs( charge(ion) ) < 0.1 )
	{
	    throw std::runtime_error( "Error: the ion you chose is in fact netural, thus add ion won't work" );
	}

	double buffer = ion_radius + max_radius + shlext;
        std::cout << "Info: rmax, rion, shell: " << max_radius << " " << ion_radius << " " << shlext << std::endl;

        std::cout << "Info: pmin: " << pmin[0] << " " << pmin[1] << " " << pmin[2] << std::endl;
	std::cout << "Info: pmax: " << pmax[0] << " " << pmax[1] << " " << pmax[2] << std::endl;


        pmin -= makevec(buffer, buffer, buffer);
        pmax += makevec(buffer, buffer, buffer);
        std::cout << "Info: enclosing: " << pmin[0] << " " << pmin[1] << " " << pmin[2] << std::endl;
	std::cout << "Info:        to: " << pmax[0] << " " << pmax[1] << " " << pmax[2] << std::endl;


        double edge = max(pmax - pmin);
        int ndepth = int( log(edge/resolution)/log(2.0) ) + 1;
        octree_t tree( pmin, ndepth, resolution );


	double outer_radius = max_radius + ion_radius + shlext;
        double closeness = 1.0;
   
        solvent_shaper shaper(ionee, closeness, outer_radius);
        tree.make_shape( shaper, ion_radius );

        std::pair<double, numvec> min = make_pair( 1e20, numvec(3));
	std::pair<double, numvec> max = make_pair(-1e20, numvec(3));
        tree.calculate( capbox::elepot_t(&pcord[0], &pchrg[0], solute_natom), min, max );

        for( int i=0; i < num_ion; ++i )
        {
            numvec pos = charge(ion)>0 ? min.second : max.second;

	    if( nosolvent )
            {
                // no solvent, simply add the solvent to the end
                std::cout << "Info: new ion will be placed at: ";
                std::cout << pos[0] << " " << pos[1] << " "  << pos[2] << " " << std::endl;

                numvec sft = pos - center(ion);
                ionee.insert_resd( ion, sft, -1 );
            }
            else
	    {
                int rid=-1;
                if( find_bump_solvent(ionee, pos, ion_radius, rid) )
                {
                    // solvent present, find the nearest solvent residue, use ion to replace it.
                    // besides the ion must be  before solvent, after the solute
                    numvec sft = pos - center(ion);
                    ionee.insert_resd( ion, sft, solute_nresd+i );
                    ionee.remove_resd( rid );
                }
                else
                {
                    // no bump atom: insert residue directly, no need to delete solvent residue
                    numvec sft = pos - center(ion);
                    ionee.insert_resd( ion, sft, solute_nresd+i );
                }
            }

            if( i < num_ion-1 )
            {
                // more ions to add: cut a sphere from the tree, make sure no overlapping
                // the radius is 2*ion_radius, to avoid in the future ion
                // with inverse sign will be placed close to this ion.
                outsphere_shaper shaper2(pos, ion_radius*2);
                tree.make_shape( shaper2, ion_radius );

                min.first = 1e20;
                max.first =-1e20;
                tree.calculate(elepot_t(pos, charge(ion)), min, max);
            }
        }

/*
        solute_size += num_ion*ion.natom();

        if( solute_size < mol.natom() )
        {
            mol.set_i(SOLUTE_NATOM, solute_size);
        }
 */       
    }

} // namespace mort

