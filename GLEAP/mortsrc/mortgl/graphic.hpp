#ifndef MORT_GUILIB_GRAPHIC_HPP
#define MORT_GUILIB_GRAPHIC_HPP

#include <object.hpp>

#include "glutils.hpp"

namespace mort
{
    class glstyle_t;

    class graphic_t
    {
    public:

        graphic_t() {}

        graphic_t(const molecule_t& m);

	virtual ~graphic_t() {}

        virtual string name() const = 0;

        virtual void paint() = 0;

	virtual void update() = 0;

        virtual shared_ptr<graphic_t> clone(const molecule_t& m, const vector<string>& args) const = 0;

        void setall(const hashid_t& parmid, const hashid_t& value);

        void set(const atomvec_t& bound, const hashid_t& parmid, const hashid_t& value);

	bool get(int absid, const hashid_t& parmid, hashid_t& v) const;

 	static shared_ptr<graphic_t> construct(const molecule_t& m, const string& type, const string& argument);

	static std::map<string, const graphic_t*>& factory();

    private:

        map< hashid_t, map<int, hashid_t> > m_params;

	const molecule_t* m_pmol;
    };

    typedef shared_ptr<graphic_t> graphic_ptr;


    class molecule_graphic : public graphic_t
    {
    public:
 
        molecule_graphic();

        molecule_graphic(const molecule_t& mol, const hashid_t& style);

        virtual ~molecule_graphic();

        virtual string name() const;

        virtual void paint();

	virtual void update();

	virtual shared_ptr<graphic_t> clone(const molecule_t& m, const vector<string>& args) const;

    private:

        const molecule_t*  m_pmol;

        string m_name;

        shared_ptr< glstyle_t > m_style;

    };

    class label_graphic : public graphic_t
    {
    public:

         label_graphic();
      
         label_graphic(const molecule_t& mol, const hashid_t& level, const hashid_t& parm);

	 virtual ~label_graphic();

	 virtual string name() const;

	 virtual void paint();

	 virtual void update();

	 virtual shared_ptr<graphic_t> clone(const molecule_t& m, const vector<string>& args) const;

    private:

        const molecule_t* m_pmol;

	hashid_t m_level;

	hashid_t m_parm;

	string m_name;

	int m_base;
    };

    class ribbon_graphic : public graphic_t
    {
    public:

        ribbon_graphic();

        ribbon_graphic(const molecule_t& pmol );

        virtual ~ribbon_graphic( );

        virtual string name() const;

        virtual void paint();

	virtual void update();

	virtual shared_ptr<graphic_t> clone(const molecule_t& m, const vector<string>& args) const;

    private:

        string m_name;

        const molecule_t* m_pmol;
    };
    
    class surface_graphic : public graphic_t
    {
    public:

        surface_graphic(const molecule_t& m);

        virtual ~surface_graphic();

        virtual string name() const;

        virtual void paint();

	virtual void update();

        virtual shared_ptr<graphic_t> clone(const molecule_t& m, const vector<string>& args) const;

    private:

        string m_name;
 
        shared_ptr< molecule_t > m_molecule;
    };

}
        
#endif

