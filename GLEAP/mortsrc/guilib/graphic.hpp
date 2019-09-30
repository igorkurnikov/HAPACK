#ifndef MORT_GUILIB_GRAPHIC_HPP
#define MORT_GUILIB_GRAPHIC_HPP

#include <object.hpp>

#include "glutils.hpp"

namespace mort
{
    class style_t;

    class graphic_i
    {
    public:
        graphic_i() {}

	virtual ~graphic_i() {}

        virtual string name() const = 0;

        virtual void paint() = 0;
    };

    class molecule_graphic : public graphic_i
    {
    public:

        molecule_graphic( const molecule_t& mol, const hashid_t& style );

        virtual ~molecule_graphic();

        virtual string name() const;

        virtual void paint();

    private:

        const molecule_t&  m_molecule;

        string m_name;

        shared_ptr< style_t > m_style;

    };

    class ribbon_graphic : public graphic_i
    {
    public: 

        ribbon_graphic( molecule_t& pmol );

        virtual ~ribbon_graphic( );

        virtual string name() const;

        virtual void paint();

    private:

        string m_name;

        molecule_t* m_pmol;
    };
    
    class surface_graphic : public graphic_i
    {
    public:

        surface_graphic( shared_ptr< molecule_t > pmol );

        virtual ~surface_graphic();

        virtual string name() const;

        virtual void paint();

    private:

        string m_name;
 
        shared_ptr< molecule_t > m_molecule;
    };

}
        
#endif

