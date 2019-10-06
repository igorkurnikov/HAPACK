#ifndef GLEAP_MORTSRC_CAPBOX_SHAPER_HPP
#define GLEAP_MORTSRC_CAPBOX_SHAPER_HPP

#include <common.hpp>

namespace mort
{
    class molecule_t;
 
    class ionee_i;
   
    namespace capbox
    {

        enum shaper_result_e { UNKNOWN_TYPE, INCLUDED, PARTIAL, EXCLUDED };

        class shaper_i
        {
        public:

            shaper_i() {}

            virtual ~shaper_i() {}

            virtual int check_atm( const numvec& pos, double rvdw ) = 0;
            
            virtual int check_box( const numvec& center, const numvec& size, double rvdw ) = 0;
        };

        typedef shared_ptr< shaper_i > shaper_ptr;

        class solvent_shaper : public shaper_i
        {
        public:

            solvent_shaper( ionee_i& solute, double closeness, double outer_radius );

            virtual ~solvent_shaper() {}

            virtual int check_atm( const numvec& pos, double rvdw );

            virtual int check_box( const numvec& center, const numvec& size, double rsolvent );

        private:

            int check_each( const numvec& center, double radius, double rsolvent, int isphere ) const;

            void make_exclusion_list( double radius );

            int m_natom;

            vector<double>& m_coord;

            vector<double>& m_pvdwr;

            double m_closeness;

            double m_outer_radius;

            vector<int> m_exclusion;

            vector<double> m_exclusion_radii;
        
        };

        class outsphere_shaper : public shaper_i
        {
        public:

            outsphere_shaper( const numvec& sphere, double radius );

            virtual ~outsphere_shaper() {}

            virtual int check_atm( const numvec& pos, double rvdw );

            virtual int check_box( const numvec& center, const numvec& size, double rsolvent );

        private:

            numvec m_sphere;

            double m_radius;
        };


    


    } // namespace capbox

} // namespace mort


#endif

