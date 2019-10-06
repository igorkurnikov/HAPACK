#ifndef MORT_ANYSTR_H
#define MORT_ANYSTR_H

#include <map>
#include <vector>
#include <string>
#include <iosfwd>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <common/numvec.hpp>
#include <common/hashcode.hpp>


namespace mort
{
    using std::map;
    using std::string;
    using std::ostream;
    using boost::any;
    using boost::shared_ptr;

    /// \brief entity_t is a container of parameters.
    /// \ingroup objdef
    /// entity_t provides two set of methods to access parameters
    /// firstly, you can access (set or get) parameters of "any" type 
    ///    by defining a parameter type.
    /// secondly, you can always access parameter of type string if the
    ///    parameter is not very common.
    class root_t
    {
    public:

        /// \brief constructor
        root_t();
        
        /// \brief copy constructor
        root_t( const root_t& rhs );
        
        /// \brief assignment opertaor
        root_t& operator=( const root_t& rhs );
        
        /// \brief deconstructor
        virtual ~root_t();
        
        /// \brief swap the content of two entity
        ///
        /// very useful for exceptional safety.
        virtual void swap( root_t& rhs );
        
        /// \brief setter function       
        void set_i(const hashid_t& parmid, const int& value); 
        void set_d(const hashid_t& parmid, const double& value); 
	void set_s(const hashid_t& parmid, const string& value); 
	void set_v(const hashid_t& parmid, const numvec& value); 
	void set_a(const hashid_t& parmid, const any& value); 
	
	void set_i(const string& parmname, const int& value);
        void set_d(const string& parmname, const double& value);
	void set_s(const string& parmname, const string& value);
	void set_v(const string& parmname, const numvec& value);
	void set_a(const string& parmname, const any& value);

        // type I getter function
        bool get_i(const hashid_t& parmid, int& value) const;
	bool get_d(const hashid_t& parmid, double& value) const;
	bool get_s(const hashid_t& parmid, string& value) const;
	bool get_v(const hashid_t& parmid, numvec& value) const;
	bool get_a(const hashid_t& parmid, any& value) const;

        bool get_i(const string& parmname, int& value) const;
	bool get_d(const string& parmname, double& value) const;
	bool get_s(const string& parmname, string& value) const;
	bool get_v(const string& parmname, numvec& value) const;
	bool get_a(const string& parmname, any& value) const;

        bool get_iptr(const hashid_t& parmid, int*& value);
	bool get_dptr(const hashid_t& parmid, double*& value);
	bool get_sptr(const hashid_t& parmid, string*& value);
	bool get_vptr(const hashid_t& parmid, numvec*& value);
	bool get_aptr(const hashid_t& parmid, any*& value);

        bool get_iptr(const string& parmname, int*& value);
	bool get_dptr(const string& parmname, double*& value);
	bool get_sptr(const string& parmname, string*& value);
	bool get_vptr(const string& parmname, numvec*& value);
	bool get_aptr(const string& parmname, any*& value);

	bool get_iptr(const hashid_t& parmid, const int*& value) const;
	bool get_dptr(const hashid_t& parmid, const double*& value) const;
	bool get_sptr(const hashid_t& parmid, const string*& value) const;
	bool get_vptr(const hashid_t& parmid, const numvec*& value) const;
	bool get_aptr(const hashid_t& parmid, const any*& value) const;

	bool get_iptr(const string& parmname, const int*& value) const;
	bool get_dptr(const string& parmname, const double*& value) const;
	bool get_sptr(const string& parmname, const string*& value) const;
	bool get_vptr(const string& parmname, const numvec*& value) const;
	bool get_aptr(const string& parmname, const any*& value) const;

        // type II getter function
        int& get_i(const hashid_t& parmid);
        double& get_d(const hashid_t& parmid);
        string& get_s(const hashid_t& parmid);
        numvec& get_v(const hashid_t& parmid);
        any& get_a(const hashid_t& parmid);

        int& get_i(const string& parmname);
        double& get_d(const string& parmname);
        string& get_s(const string& parmname);
        numvec& get_v(const string& parmname);
        any& get_a(const string& parmname);

        const int& get_i(const hashid_t& parmid) const;
        const double& get_d(const hashid_t& parmid) const;
        const string& get_s(const hashid_t& parmid) const;
        const numvec& get_v(const hashid_t& parmid) const;
        const any& get_a(const hashid_t& parmid) const;

        const int& get_i(const string& parmname) const;
        const double& get_d(const string& parmname) const;
        const string& get_s(const string& parmname) const;
        const numvec& get_v(const string& parmname) const;
        const any& get_a(const string& parmname) const;

        /// type III getter function
	int& frcget_i(const hashid_t& parmid); 
        double& frcget_d(const hashid_t& parmid);
        string& frcget_s(const hashid_t& parmid);
        numvec& frcget_v(const hashid_t& parmid);
        any& frcget_a(const hashid_t& parmid);

        int& frcget_i(const string& parmname); 
        double& frcget_d(const string& parmname);
        string& frcget_s(const string& parmname);
        numvec& frcget_v(const string& parmname);
        any& frcget_a(const string& parmname);


        // Testing functions
        bool has_i( const hashid_t& parmid ) const;
	bool has_d( const hashid_t& parmid ) const;
	bool has_s( const hashid_t& parmid ) const;
	bool has_v( const hashid_t& parmid ) const;
	bool has_a( const hashid_t& parmid ) const;

        bool has_i( const string& parmid ) const;
	bool has_d( const string& parmid ) const;
	bool has_s( const string& parmid ) const;
	bool has_v( const string& parmid ) const;
	bool has_a( const string& parmid ) const;

        // Iterator functions
        map<hashid_t, int>::const_iterator ibegin() const;
        map<hashid_t, int>::const_iterator iend() const;

	map<hashid_t, double>::const_iterator dbegin() const;
	map<hashid_t, double>::const_iterator dend() const;

	map<hashid_t, string>::const_iterator sbegin() const;
	map<hashid_t, string>::const_iterator send() const;

	map<hashid_t, numvec>::const_iterator vbegin() const;
	map<hashid_t, numvec>::const_iterator vend() const;

	map<hashid_t, any>::const_iterator abegin() const;
	map<hashid_t, any>::const_iterator aend() const;

      private:

        map<hashid_t, int> m_icontent;
	map<hashid_t, double> m_dcontent;
	map<hashid_t, string> m_scontent;
	map<hashid_t, numvec> m_vcontent;
	map<hashid_t, any> m_acontent;
    };    
     
    ostream& mortlog( int level );

    typedef shared_ptr<root_t> entity_ptr;

    bool ismol( const entity_ptr& pe );

    bool ismdb( const entity_ptr& pe );
    
    bool isavec( const entity_ptr& pe );
    
    bool isbvec( const entity_ptr& pe );

} // namespace mort

#endif

