#ifndef MORTSRC_OBJECT_OBJDEF_HPP
#define MORTSRC_OBJECT_OBJDEF_HPP

#include <string>
#include <boost/any.hpp>
#include <common.hpp>
#include "mfwd.hpp"

namespace mort
{
    using std::string;
    using boost::any;

    /// \brief morf_t provide a convinent way to access parameter
    /// \ingroup objdef
    /// To access parameter through component_t is inconvinent,
    /// morf_t packs the related function together and make it much
    /// easier. A lot of types are derived from it, such as atom_i, 
    /// bond_i, etc.
    class morf_t
    {
    public:

        /// \brief constructor
        morf_t( const molecule_t& mol, int compid, int objectid );
        morf_t( const morf_t& rhs );
        virtual ~morf_t();
        morf_t& operator=( const morf_t& rhs );
   

        /// ID functions
        int absid() const;  /// get absolute ID
        int relid() const;  /// get relative ID
        hashid_t cmpid() const;  /// get component ID

        /// \brief return the reference of the molecule it belongs to
        molecule_t& getmol();
        molecule_t const& getmol() const;
        mcmpdata_t* getcmp(); 
        mcmpdata_t const* getcmp() const;
        mcmprela_t* getadj(int cid); 
	mcmprela_t const* getadj(int cid) const;

        // setter function
	void set_i( const hashid_t& parmid, const int& value ); 
        void set_d( const hashid_t& parmid, const double& value );
        void set_s( const hashid_t& parmid, const string& value );
        void set_v( const hashid_t& parmid, const numvec& value );
        void set_a( const hashid_t& parmid, const any& value );

        void set_i( const string& parmname, const int& value );
	void set_d( const string& parmname, const double& value ); 
        void set_s( const string& parmname, const string& value ); 
	void set_v( const string& parmname, const numvec& value );
	void set_a( const string& parmname, const any& value ); 

        // type I getter function
	bool get_i( const hashid_t& parmid, int& value ) const; 
	bool get_d( const hashid_t& parmid, double& value ) const;
	bool get_s( const hashid_t& parmid, string& value ) const;
	bool get_v( const hashid_t& parmid, numvec& value ) const;
	bool get_a( const hashid_t& parmid, any& value ) const; 

	bool get_i( const string& parmname, int& value ) const; 
	bool get_d( const string& parmname, double& value ) const;
	bool get_s( const string& parmname, string& value ) const;
	bool get_v( const string& parmname, numvec& value ) const;
	bool get_a( const string& parmname, any& value ) const;

        bool get_iptr(const hashid_t& parmid, int*& pi);
        bool get_dptr(const hashid_t& parmid, double*& pd);
        bool get_sptr(const hashid_t& parmid, string*& ps);
        bool get_vptr(const hashid_t& parmid, double*& pv);
        bool get_aptr(const hashid_t& parmid, any*& pa);

        bool get_iptr(const hashid_t& parmid, const int*& pi) const;
	bool get_dptr(const hashid_t& parmid, const double*& pd) const;
        bool get_sptr(const hashid_t& parmid, const string*& ps) const;
        bool get_vptr(const hashid_t& parmid, const double*& pv) const;
        bool get_aptr(const hashid_t& parmid, const any*& pa) const;

	/// type II getter function
        int get_i(const hashid_t& parmid) const;
        double get_d(const hashid_t& parmid) const;
        string get_s(const hashid_t& parmid) const;
	numvec get_v(const hashid_t& parmid) const;
        any get_a(const hashid_t& parmid) const;

        int get_i(const string& parmname) const;
        double get_d(const string& parmname) const;
        string get_s(const string& parmname) const;
	numvec get_v(const string& parmname) const;
        any get_a(const string& parmname) const;

        int* get_iptr(const hashid_t& parmid);
        double* get_dptr(const hashid_t& parmid);
        string* get_sptr(const hashid_t& parmid);
        double* get_vptr(const hashid_t& parmid);
        any* get_aptr(const hashid_t& parmid);

        const int* get_iptr(const hashid_t& parmid) const;
        const double* get_dptr(const hashid_t& parmid) const;
        const string* get_sptr(const hashid_t& parmid) const;
	const double* get_vptr(const hashid_t& parmid) const;
        const any* get_aptr(const hashid_t& parmid) const;

        /// Testing function
        bool has_i(const hashid_t& parmid) const;
        bool has_d(const hashid_t& parmid) const;
        bool has_s(const hashid_t& parmid) const;
        bool has_v(const hashid_t& parmid) const;
        bool has_a(const hashid_t& parmid) const;

        bool has_i(const string& parmname) const; 
        bool has_d(const string& parmname) const; 
        bool has_s(const string& parmname) const; 
	bool has_v(const string& parmname) const; 
        bool has_a(const string& parmname) const; 

        bool has_iptr(const hashid_t& parmid) const;
        bool has_dptr(const hashid_t& parmid) const;
        bool has_sptr(const hashid_t& parmid) const;
        bool has_vptr(const hashid_t& parmid) const;
        bool has_aptr(const hashid_t& parmid) const;


        // connection functions
	bool connect(const morf_t& ptr);
	bool is_connected_to( const morf_t& ptr) const;
        bool disconnect(const morf_t& ptr);

        // iterator functions
        mobjiter_t mobj_begin( const hashid_t& id ) const;
        mobjiter_t mobj_end  ( const hashid_t& id ) const;
        atomiter_t atom_begin() const;
        atomiter_t atom_end  () const;
        bonditer_t bond_begin() const;
        bonditer_t bond_end  () const;
        resditer_t resd_begin() const;
	resditer_t resd_end  () const;
        angliter_t angl_begin() const;
        angliter_t angl_end  () const;
        diheiter_t dihe_begin() const;
        diheiter_t dihe_end  () const;
        impriter_t impr_begin() const;
	impriter_t impr_end  () const;
        mobjiter_t tor2_begin() const;
        mobjiter_t tor2_end  () const;
        mobjiter_t ptor_begin() const;
        mobjiter_t ptor_end  () const;
        

        std::vector<int> const& related_atom_ids() const;
        std::vector<int> const& related_bond_ids() const;
        std::vector<int> const& related_resd_ids() const;

        // range functions
        mobj_range mobjs(const hashid_t& id) const;
        atom_range atoms() const; 
        bond_range bonds() const;
        angl_range angls() const;
        dihe_range dihes() const;
        resd_range resds() const;
        mobj_range tor2s() const; 
        mobj_range ptors() const;
 
        int nmobj(const hashid_t& id) const;
        int natom() const { return nmobj(ATOM); }
	int nbond() const { return nmobj(BOND); }
	int nresd() const { return nmobj(RESD); }
	int nangl() const { return nmobj(ANGL); }
	int noops() const { return nmobj(OOPS); }
	int ntors() const { return nmobj(TORS); }
	int ntor2() const { return nmobj(TOR2); }
	int nptor() const { return nmobj(PTOR); }

	template<typename T> friend class iter_T;

    private:

        molecule_t* m_molecule;
        hashid_t m_cmpid;
        int m_absid;
    };

    /// \brief operator ==
    ///
    /// compare the ID number of two object; 
    bool operator ==( const morf_t& lhs, const morf_t& rhs );
 
    /// \brief operator !=
    ///
    /// compare the ID number of two object;
    bool operator != ( const morf_t& lhs, const morf_t& rhs );

    /// \brief operator <
    ///
    /// compare the ID number of two object; 
    bool operator < ( const morf_t& lhs, const morf_t& rhs );

    /// \brief operator <
    ///
    /// compare the ID number of two object; 
    bool operator > ( const morf_t& lhs, const morf_t& rhs ); 

    /// \brief operator <
    ///
    /// compare the ID number of two object; 
    bool operator <=( const morf_t& lhs, const morf_t& rhs );

    /// \brief operator <
    ///
    /// compare the ID number of two object; 
    bool operator >=( const morf_t& lhs, const morf_t& rhs ); 

} // namespace mort


#endif

