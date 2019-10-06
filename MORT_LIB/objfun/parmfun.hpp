#ifndef MORT_OBJECT_PARMFUN_HPP
#define MORT_OBJECT_PARMFUN_HPP

#include <iosfwd>
#include <string>
#include <vector>
#include <boost/any.hpp>
#include <common/numvec.hpp>
#include <common/hashcode.hpp>

namespace mort
{
    using std::vector;
    using std::string;
    using boost::any;

    class morf_t;

    class molecule_t;

    struct iparm_setter2
    {
        typedef void result_type;
        iparm_setter2(const hashid_t& parmid);
	void operator()(morf_t& p, int value) const;
	hashid_t m_parmid;
    };

    struct dparm_setter2
    {
        typedef void result_type;
        dparm_setter2(const hashid_t& parmid);
	void operator()(morf_t& p, double value) const;
	hashid_t m_parmid;
    };

    struct sparm_setter2
    {
        typedef void result_type;
        sparm_setter2(const hashid_t& parmid);
	void operator()(morf_t& p, const string& value) const;
	hashid_t m_parmid;
    };

    struct vparm_setter2
    {
        typedef void result_type;
        vparm_setter2(const hashid_t& parmid);
	void operator()(morf_t& p, const numvec& value) const;
	hashid_t m_parmid;
    };

    struct aparm_setter2
    {
        typedef void result_type;
        aparm_setter2(const hashid_t& parmid);
	void operator()(morf_t& p, const any& value) const;
	hashid_t m_parmid;
    };

    struct iparm_setter1
    {
        typedef void result_type;
        iparm_setter1(const hashid_t& parmid, int value);
	void operator()(morf_t& p) const;
	hashid_t m_parmid;
	int m_value;
    };

    struct dparm_setter1
    {
        typedef void result_type;
        dparm_setter1(const hashid_t& parmid, double value);
	void operator()(morf_t& p) const;
	hashid_t m_parmid;
	double m_value;
    };

    struct sparm_setter1
    {
        typedef void result_type;
        sparm_setter1(const hashid_t& parmid, const string& value);
	void operator()(morf_t& p) const;
	hashid_t m_parmid;
	string m_value;
    };

    struct vparm_setter1
    {
        typedef void result_type;
        vparm_setter1(const hashid_t& parmid, const numvec& value);
	void operator()(morf_t& p) const;
	hashid_t m_parmid;
	numvec m_value;
    };

    struct aparm_setter1
    {
        typedef void result_type;
        aparm_setter1(const hashid_t& parmid, const any& value);
	void operator()(morf_t& p) const;
	hashid_t m_parmid;
	any m_value;
    };

    struct iparm_getter
    {
        iparm_getter(const hashid_t& parmid);
	int operator()(const morf_t& p) const;
	hashid_t m_parmid;
    };

    struct dparm_getter
    {
        dparm_getter(const hashid_t& parmid);
	double operator()(const morf_t& p) const;
	hashid_t m_parmid;
    };

    struct sparm_getter
    {
        sparm_getter(const hashid_t& parmid);
	string operator()(const morf_t& p) const;
	hashid_t m_parmid;
    };

    struct vparm_getter
    {
        vparm_getter(const hashid_t& parmid);
	numvec operator()(const morf_t& p) const;
	hashid_t m_parmid;
    };

    struct aparm_getter
    {
        aparm_getter(const hashid_t& parmid);
	any operator()(const morf_t& p) const;
	hashid_t m_parmid;
    };

    struct iparm_cmper1
    {
        typedef bool result_type;
	typedef const morf_t& argument_type;
	iparm_cmper1(const hashid_t& parmid, int value);
	bool operator()(const morf_t& p) const;

        hashid_t m_parmid;
	int m_value;
    };

    struct dparm_cmper1
    {
        typedef bool result_type;
	typedef const morf_t& argument_type;
	dparm_cmper1(const hashid_t& parmid, double value);
	bool operator()(const morf_t& p) const;

        hashid_t m_parmid;
	double m_value;
    };

    struct sparm_cmper1
    {
        typedef bool result_type;
	typedef const morf_t& argument_type;
	sparm_cmper1(const hashid_t& parmid, const string& value);
	bool operator()(const morf_t& p) const;

        hashid_t m_parmid;
	string m_value;
    };

    struct vparm_cmper1
    {
        typedef bool result_type;
	typedef const morf_t& argument_type;
	vparm_cmper1(const hashid_t& parmid, const numvec& value);
	bool operator()(const morf_t& p) const;

        hashid_t m_parmid;
	numvec m_value;
    };

    struct iparm_cmper2
    {
        typedef bool result_type;
	typedef const morf_t& first_argument_type;
	typedef int second_argument_type;

        iparm_cmper2(const hashid_t& parmid);
	bool operator()(const morf_t& p, int value) const;

	hashid_t m_parmid;
    };

    struct dparm_cmper2
    {
        typedef bool result_type;
	typedef const morf_t& first_argument_type;
	typedef int second_argument_type;

        dparm_cmper2(const hashid_t& parmid);
	bool operator()(const morf_t& p, double value) const;

	hashid_t m_parmid;
    };

    struct sparm_cmper2
    {
        typedef bool result_type;
	typedef const morf_t& first_argument_type;
	typedef int second_argument_type;

        sparm_cmper2(const hashid_t& parmid);
	bool operator()(const morf_t& p, const string& value) const;

	hashid_t m_parmid;
    };

    struct vparm_cmper2
    {
        typedef bool result_type;
	typedef const morf_t& first_argument_type;
	typedef int second_argument_type;

        vparm_cmper2(const hashid_t& parmid);
	bool operator()(const morf_t& p, const numvec& value) const;

	hashid_t m_parmid;
    };

    using std::istream;
    using std::ostream;

    void read_iparm(istream& is, morf_t& mo, const hashid_t& parmid);
    void read_dparm(istream& is, morf_t& mo, const hashid_t& parmid);
    void read_sparm(istream& is, morf_t& mo, const hashid_t& parmid);
    void read_vparm(istream& is, morf_t& mo, const hashid_t& parmid, int nitem);

    void write_iparm(ostream& os, const morf_t& mo, const hashid_t& parmid, const string& format);
    void write_dparm(ostream& os, const morf_t& mo, const hashid_t& parmid, const string& format);
    void write_sparm(ostream& os, const morf_t& mo, const hashid_t& parmid, const string& format);
    void write_vparm(ostream& os, const morf_t& mo, const hashid_t& parmid, const string& format);

    class fortran_t;
    void fortran_write_iparm( fortran_t* format, morf_t& p, const hashid_t& parmid );
    void fortran_write_dparm( fortran_t* format, morf_t& p, const hashid_t& parmid );
    void fortran_write_sparm( fortran_t* format, morf_t& p, const hashid_t& parmid );


    class root_t;
    bool copy_iparm( const morf_t& src, const hashid_t& parmid, morf_t& dst );
    bool copy_dparm( const morf_t& src, const hashid_t& parmid, morf_t& dst );
    bool copy_sparm( const morf_t& src, const hashid_t& parmid, morf_t& dst );
    bool copy_vparm( const morf_t& src, const hashid_t& parmid, morf_t& dst );
    bool copy_aparm( const morf_t& src, const hashid_t& parmid, morf_t& dst );
    void copy_allparms(const morf_t& src, morf_t& dst);
    void copy_allparms_except(const morf_t& src, morf_t& dst, const vector<hashid_t>& excepts);   
    void copy_allparms_except(const morf_t& src, morf_t& dst, const hashid_t& except);   
    void copy_allparms( const root_t& src, root_t& dst );

    string parm2str( const morf_t& mo, const hashid_t& parmid );

    vector<int>&    get_ivec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid);
    vector<double>& get_dvec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid);
    vector<string>& get_svec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid);
    vector<double>& get_vvec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid);

    vector<int>    const& get_ivec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid);
    vector<double> const& get_dvec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid);
    vector<string> const& get_svec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid);
    vector<double> const& get_vvec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid);

    double charge( const molecule_t& m );

    double charge( const morf_t& mo );

    double weight( const molecule_t& m );

    double weight( const morf_t& mo );

    double get_vdwr( const string& type );

    void set_vdwr( const molecule_t& m, vector<double>& vdwr );

    void set_vdwr( molecule_t& m );

    void set_vdwr( morf_t& mo );

} // namespace mort

#endif

