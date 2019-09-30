#ifndef DDUTILS_MORT_COMPONENT_H
#define DDUTILS_MORT_COMPONENT_H

#include <set>
#include <map>
#include <vector>
#include <string>
#include <boost/any.hpp>
#include <common.hpp>

namespace mort
{
    using std::set;
    using std::map;
    using std::vector;
    using std::string; 
    using boost::any;

    typedef map<hashid_t, vector<int> >::const_iterator iparm_iterator;
    typedef map<hashid_t, vector<double> >::const_iterator dparm_iterator;
    typedef map<hashid_t, vector<string> >::const_iterator sparm_iterator;
    typedef map<hashid_t, vector<double> >::const_iterator vparm_iterator;
    typedef map<hashid_t, vector<any> >::const_iterator aparm_iterator;

    /// \brief component_t contain the parameters of all the objects
    ///
    /// component_t is included by molecule_t to hold all the parameters.
    /// it provides the inteface to set parameter of each individual object
    /// the relationships between components are recorded in class adjacency
    ///
    /// \sa adjacency_t, molecule_t
    class mcmpdata_t
    {
      public:

        /// \briefthe iterator type of component,
        ///
        /// which is acturally of type vector< int >::iterator
        /// Dereferencing it will give you an integer which is ID number of an object.
        typedef vector<int>::iterator iterator;

        /// \brief contructor
        /// \param size the number of objects in the component;
        mcmpdata_t(int size = 0);
        
        /// \brief copy constructor
        mcmpdata_t(const mcmpdata_t& rhs);

        /// \brief deconstructor
        ~mcmpdata_t();

	/// \brief swap the content of two component
	void swap( mcmpdata_t& rhs );
        
        /// \brief add new object to the end.
        iterator append(int n=1);

        /// \brief insert new object at the certain position.
        /// \param i the position where new object will be added.
        iterator insert(int relid, int n=1);

        /// \brief remove object with given ID number
        /// \param the ID number of the object going to be removed.
        bool remove(int id);
        
        /// \brief adjust the size of the component
        /// \param size the new size
        void resize(int size);

        /// \brief return the size of the component
        int size() const;

        iterator absid_begin() const; /// \brief the begin point of objects
        iterator absid_end() const; /// \brief the ending point of objects

        bool isclean() const; // test if the component is clean, i.e. contineous absolute ID
        void cleanup(map<int, int>& old2new); // clean up the component, return an old to new ID map.

        // setter function
        void set_i(const hashid_t& parmid, int absid, const int& value);  
        void set_d(const hashid_t& parmid, int absid, const double& value); 
        void set_s(const hashid_t& parmid, int absid, const string& value); 
        void set_v(const hashid_t& parmid, int absid, const numvec& value); 
        void set_a(const hashid_t& parmid, int absid, const any& value);

        // type I getter function
        bool get_i(const hashid_t& parmid, int absid, int& v) const; 
        bool get_d(const hashid_t& parmid, int absid, double& v) const;
	bool get_s(const hashid_t& parmid, int absid, string& v) const;
        bool get_v(const hashid_t& parmid, int absid, numvec& v) const;
        bool get_a(const hashid_t& parmid, int absid, any& v) const;

        // get pointer to parameter
        bool get_iptr(const hashid_t& parmid, int absid, int*& pi);
        bool get_dptr(const hashid_t& parmid, int absid, double*& pd);
        bool get_sptr(const hashid_t& parmid, int absid, string*& ps);
        bool get_vptr(const hashid_t& parmid, int absid, double*& pv);
        bool get_aptr(const hashid_t& parmid, int absid, any*& pa);

        bool get_iptr(const hashid_t& parmid, int absid, const int*& pi) const;
	bool get_dptr(const hashid_t& parmid, int absid, const double*& pd) const;
        bool get_sptr(const hashid_t& parmid, int absid, const string*& ps) const;
        bool get_vptr(const hashid_t& parmid, int absid, const double*& pv) const;
        bool get_aptr(const hashid_t& parmid, int absid, const any*& pa) const;

        // type II getter function
        int& get_i( const hashid_t& parmid, int absid );
        double& get_d( const hashid_t& parmid, int absid );
        string& get_s( const hashid_t& parmid, int absid );
        any& get_a( const hashid_t& parmid, int absid );

        const int& get_i( const hashid_t& parmid, int absid ) const;
        const double& get_d( const hashid_t& parmid, int absid ) const;
        const string& get_s( const hashid_t& parmid, int absid ) const;
        numvec get_v( const hashid_t& parmid, int absid ) const;
        const any& get_a( const hashid_t& parmid, int absid ) const;

        // get pointer to parameter
        double* get_dptr( const hashid_t& parmid, int absid );
        double* get_vptr( const hashid_t& parmid, int absid );
	string* get_sptr( const hashid_t& parmid, int absid );

	double const* get_dptr(const hashid_t& parmid, int absid) const;
	double const* get_vptr(const hashid_t& parmid, int absid) const; 
        string const* get_sptr(const hashid_t& parmid, int absid) const;

        // get reference to holder
        vector<int>& get_ivec(const hashid_t& parmid);
        vector<double>& get_dvec(const hashid_t& parmid);
        vector<double>& get_vvec(const hashid_t& parmid);
        vector<string>& get_svec(const hashid_t& parmid);

        vector<int> const& get_ivec(const hashid_t& parmid) const;
        vector<double> const& get_dvec(const hashid_t& parmid) const;
        vector<double> const& get_vvec(const hashid_t& parmid) const;
        vector<string> const& get_svec(const hashid_t& parmid) const;

        // Testing functions
        bool has_i( const hashid_t& parmid, int absid ) const; 
        bool has_d( const hashid_t& parmid, int absid ) const; 
        bool has_s( const hashid_t& parmid, int absid ) const; 
        bool has_v( const hashid_t& parmid, int absid ) const; 
        bool has_a( const hashid_t& parmid, int absid ) const; 

        // Iterator functions
        iparm_iterator iparm_begin() const;
        iparm_iterator iparm_end()   const;
        dparm_iterator dparm_begin() const;
	dparm_iterator dparm_end()   const;
        sparm_iterator sparm_begin() const;
        sparm_iterator sparm_end()   const;
        vparm_iterator vparm_begin() const;
        vparm_iterator vparm_end()   const;
        aparm_iterator aparm_begin() const; 
        aparm_iterator aparm_end()   const;
 
      private:

        void resize_absids(int size);
        void resize_parm( );

      private:

        map< hashid_t, set<int> > m_ioccupy;
        map< hashid_t, set<int> > m_doccupy;
        map< hashid_t, set<int> > m_soccupy;
        map< hashid_t, set<int> > m_voccupy;
        map< hashid_t, set<int> > m_aoccupy;

        map< hashid_t, vector<int   > > m_icontent;
	map< hashid_t, vector<double> > m_dcontent;
	map< hashid_t, vector<string> > m_scontent;
	map< hashid_t, vector<double> > m_vcontent;        
	map< hashid_t, vector<any   > > m_acontent;

        map< hashid_t, int > m_vlength;

        vector<int> m_absids;

        int m_idcounter;

        bool m_clean;
    };


} // namespace mort

#endif


