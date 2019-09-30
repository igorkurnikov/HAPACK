#include <stdexcept>
#include <common.hpp>
#include <object.hpp>
#include "solute.hpp"
#include "region.hpp"
#include "solvate.hpp"
#include "solvent.hpp"

namespace mort
{

    using std::cout;

    // go over each residue of the solvent to check if it should be added or not.
    void addpart_solvent(solute_i& solute, const solvent_t& solvent, const numvec& shift, const region_i& region )
    {
        funstack_t::push( "addpart_solvent" );

        const vector<double>& cord = solvent.getcord();
        const vector<double>& vdwr = solvent.getvdwr();

        int nresd = solvent.nresd();
        int rsize = solvent.natom()/nresd;

        for( int ir=0; ir < nresd; ++ir )
        {
            bool included = true;
                
            int ia = ir*rsize;
            for( int j=0; j < rsize; ++j, ++ia )
            {  
                numvec curt(3);
                curt[0] = shift[0] + cord[ia*3  ];
                curt[1] = shift[1] + cord[ia*3+1];
                curt[2] = shift[2] + cord[ia*3+2];

                if( region.check_sphere(curt, vdwr[ia]) != INCLUDE )
                {
                    included = false;
                    break;
                }
            }

            if( included )
            {
                // -1 means add residue to the end.
                solute.insert_resd( solvent.getresd(ir), shift, -1 );
            }
        } 

        funstack_t::pop();
    }

    // add the full solvent directly.
    void addfull_solvent(solute_i& solute, const solvent_t& solvent, const numvec& shift)
    {
        funstack_t::push( "addfull_solvent" );
        solute.insert_resd( solvent.getfull(), shift, -1 );
        funstack_t::pop();
    }

    void solvate(solute_i& solute, const molecule_t& m, const region_i& region)
    {
        funstack_t::push( "solvate_core" );
        solvent_t solvent(m);
        numvec svt_center = solvent.center();
        numvec svt_extent = solvent.extent();
        double svt_maxrad = *std::max_element( solvent.getvdwr().begin(), solvent.getvdwr().end() );            

        numvec rgn_center = region.center();
        numvec rgn_extent = region.extent();

        assert( svt_extent.size()==3u );
        cout << format( "rgn size: %9.3f %9.3f %9.3f" ) % rgn_extent[0] % rgn_extent[1] % rgn_extent[2] << std::endl;
        cout << format( "svt size: %9.3f %9.3f %9.3f" ) % svt_extent[0] % svt_extent[1] % svt_extent[2] << std::endl;

        int nx = int(rgn_extent[0] / svt_extent[0]) + 1;
        int ny = int(rgn_extent[1] / svt_extent[1]) + 1;
        int nz = int(rgn_extent[2] / svt_extent[2]) + 1;
        std::cout << "nx,ny,nz: " << nx << " " << ny << " " << nz << std::endl;


        numvec rgn_corner(3);
        rgn_corner[0] = rgn_center[0] + 0.5*svt_extent[0]*nx;
        rgn_corner[1] = rgn_center[1] + 0.5*svt_extent[1]*ny;
        rgn_corner[2] = rgn_center[2] + 0.5*svt_extent[2]*nz;

        int ix = 0;
        int iy = 0;
        int iz = 0;
        int nc = nx*ny*nz;
        for( int i=0; i < nc; i++ )
        {
            numvec cur_center(3);
            cur_center[0] = rgn_corner[0] - (ix+0.5)*svt_extent[0];
            cur_center[1] = rgn_corner[1] - (iy+0.5)*svt_extent[1];
            cur_center[2] = rgn_corner[2] - (iz+0.5)*svt_extent[2];

            assert( svt_extent.size()==3u );
/*
            int type = region.check_rndbox(cur_center, svt_extent, svt_maxrad);
            if( type == INCLUDE ) 
            {
                addfull_solvent( solute, solvent, cur_center - svt_center );
            }
            else if( type == PARTIAL )
            {
                addpart_solvent( solute, solvent, cur_center - svt_center, region );
            }
*/
                addpart_solvent( solute, solvent, cur_center - svt_center, region );
 
            iz++;
            if(iz==nz) {iz=0; iy++;}
            if(iy==ny) {iy=0; ix++;}
        }

        funstack_t::pop();
    }

    numvec solvatecap_core( solute_i& solute, const molecule_t& m, const numvec& capcnt, double caprad, double closeness )
    {
        cap_region cap( capcnt, caprad );
        out_solute inn( solute, closeness );
        solvate( solute, m, and_region(cap, inn) );
        return cap.extent();
    }
    
    numvec solvatebox_core( solute_i& solute, const molecule_t& m, double buffer, double closeness )
    {
        funstack_t::push( "solvatebox" );

        numvec extent = regionlize( solute, true );
        extent += scalar_numvec(extent.size(), 2*buffer);

        box_region box( extent );
        out_solute inn( solute, closeness );
        solvate( solute, m, and_region(box, inn) );

        regionlize( solute, true );

        funstack_t::pop();
        return extent; 
    }
    
    // rotate to the best orientation, that longest distance is along x,y,z axis
    double solvateoct_core( solute_i& solute, const molecule_t& m, double buffer, double closeness )
    {
        centralize( solute );
        rotatelong( solute );

        numvec extent = regionlize( solute, true );        
	double extmax = max(extent) + 2*octbufsize(solute, buffer);
	extent = makevec(extmax, extmax, extmax);

        oct_region oct( extent );
        out_solute inn( solute, closeness );

        solvate( solute, m, and_region(oct,inn) );

	regionlize( solute, true );
        rotatewald( solute );
        return extmax; 
    }

    numvec solvateshl_core( solute_i& solute, const molecule_t& m, double shlext, double closeness )
    {
        regionlize( solute, true );
        shl_solute shl( solute, closeness, shlext );
        solvate( solute, m, shl );
        return shl.extent();
    } 

    void solvatecap(molecule_t& mol, const molecule_t& svt, const numvec& capcnt, double caprad, double closeness)
    {
        mort_solute solute(mol);
        solvatecap_core( solute, svt, capcnt, caprad, closeness );

        mol.set_i( SOLUTE, CAP );
        mol.set_v( CAP, makevec(capcnt[0],capcnt[1],capcnt[2],caprad) );
    }

    void solvatebox(molecule_t& mol, const molecule_t& svt, double buffer, double closeness)
    {
        mort_solute solute(mol);
        numvec extent = solvatebox_core( solute, svt, buffer, closeness );
      
        mol.set_i( SOLUTE, BOX );
        mol.set_v( BOX, makevec(extent[0], extent[1], extent[2], 90.0) );
    }     

    void solvateoct( molecule_t& mol, const molecule_t& svt, double buffer, double closeness )
    {
        mort_solute solute( mol );
        double extent = solvateoct_core( solute, svt, buffer, closeness );
        double octlen = extent*sqrt(0.75);

        mol.set_i( SOLUTE, BOX );
        mol.set_v( BOX, makevec(octlen, octlen, octlen, 109.471219) );
    }

    void solvateshl( molecule_t& mol, const molecule_t& svt, double shlext, double closeness )
    {
        mort_solute slt( mol );
        solvateshl_core( slt, svt, shlext, closeness );
    }

} // namespace mort

