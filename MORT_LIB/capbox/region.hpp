#ifndef MORTSRC_CAPBOX_REGION_HPP
#define MORTSRC_CAPBOX_REGION_HPP

#include <common.hpp>
#include "solute.hpp"

namespace mort
{
     enum region_relation_e {INCLUDE, OUTSIDE, PARTIAL};

     class region_i
     {
     public:

         virtual ~region_i() {}

         virtual numvec center() const = 0;

         virtual numvec extent() const = 0;

         /// check the intersection of the region with a sphere
         /// params:
         ///        sphcnt: center of the sphere
         ///        sphrad: radius of the sphere
         /// return value:
         ///        INCLUDE: the sphere is inside the region
         ///        OUTSIDE: the sphere is out of the region
         virtual int check_sphere( const numvec& sphcnt, double sphrad ) const = 0;

         /// check the intersection of the region with a rounded box
         /// params:
         ///        boxcnt: center of the sphere
         ///        extent: size of box
         ///        sphrad: radius of the round part
         /// return value:
         ///        INCLUDE: the box is inside the region
         ///        OUTSIDE: the box is out of the region
         ///        PARTIAL: the boxis partial in partial out
         virtual int check_rndbox( const numvec& sphcnt, const numvec& extent, double sphrad ) const = 0;
    };

    class box_region : public region_i
    {
    public:

        box_region( const numvec& boxcnt, const numvec& boxext )
            : m_cnt(boxcnt), m_ext(boxext), m_min(3), m_max(3)
        {
            m_min = boxcnt - 0.5*boxext;
            m_max = boxcnt + 0.5*boxext;
        }

        // if no center is given, consider it is on origin
        box_region( const numvec& boxext )
            : m_cnt(zero_numvec(3)), m_ext(boxext), m_min(3), m_max(3)
        {
            m_min = m_cnt - 0.5*boxext;
            m_max = m_cnt + 0.5*boxext;
        }

        virtual ~box_region() {}

        virtual int check_sphere( const numvec& sphcnt, double sphrad ) const;

        virtual int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const;

        virtual numvec center() const { return m_cnt; }

        virtual numvec extent() const { return m_ext; }

    private:

        numvec m_cnt;
        numvec m_ext;
        numvec m_min;
        numvec m_max;
    };

    class cap_region : public region_i
    {
    public:

        cap_region( const numvec& capcnt, double caprad )
           : m_cnt(capcnt),
             m_rad(caprad),
             m_ext(scalar_numvec(3, 2*caprad))
        {
        }

        virtual ~cap_region() {}

        virtual int check_sphere( const numvec& sphcnt, double sphrad ) const;

        virtual int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const;

        virtual numvec center() const { return m_cnt; }

        virtual numvec extent() const { return m_ext; }

    private:

        numvec m_cnt;
        double m_rad;
        numvec m_ext;
    };

    class oct_region : public region_i
    {
    public:

        oct_region( const numvec& octcnt, const numvec& octext )
            : m_cnt(octcnt), m_ext(octext)
        {
        }

        oct_region( const numvec& octext )
            : m_cnt(zero_numvec(3)), m_ext(octext)
        {
        }


        virtual ~oct_region() {}

        virtual int check_sphere( const numvec& sphcnt, double sphrad ) const;

        virtual int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const;

        virtual numvec center() const { return m_cnt; }

        virtual numvec extent() const { return m_ext; }

    private:

        numvec m_cnt;
        numvec m_ext;
    };

    class out_solute : public region_i
    {
    public:

        out_solute( solute_i& solute, double clsnss )
            : m_nsph(solute.natom()),
              m_cnts(solute.getcord()),
              m_rads(solute.getvdwr())
        {
            m_clsnss = clsnss;
        }

        virtual ~out_solute() {}

        virtual int check_sphere( const numvec& sphcnt, double sphrad ) const;

        virtual int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const;

        virtual numvec center() const { return zero_numvec(3); }

        virtual numvec extent() const { return scalar_numvec(3, DOUBLEMAX); }

        int check_sphere( const numvec& sphcnt, double sphrad, double& mind2 ) const;

        int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad, double& mind2 ) const;

    private:

        int m_nsph;
        const vector<double>& m_cnts;
        const vector<double>& m_rads;
        double m_clsnss;
    };


    class shl_solute : public region_i
    {
    public:

        shl_solute( solute_i& solute, double clsnss, double shlext );
   
        virtual ~shl_solute() {}

        virtual int check_sphere( const numvec& sphcnt, double sphrad ) const;

        virtual int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const;

        virtual numvec center() const { return m_center; }

        virtual numvec extent() const { return m_extent; }

    private:

        out_solute m_inner;
        numvec m_center;
        numvec m_extent;
        double m_shlext;
        double m_maxrad;
        double m_clsnss;
    };

    class and_region : public region_i
    {
    public:
  
        and_region( const region_i& lhs, const region_i& rhs );
 
        virtual ~and_region() {}

        virtual int check_sphere( const numvec& sphcnt, double sphrad ) const;

        virtual int check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const;

        virtual numvec center() const { return m_cnt; }

        virtual numvec extent() const { return m_ext; }

    private:

        const region_i* m_lhs;
        const region_i* m_rhs;
        numvec m_cnt;
        numvec m_ext;

    };

} // namespace mort

#endif

