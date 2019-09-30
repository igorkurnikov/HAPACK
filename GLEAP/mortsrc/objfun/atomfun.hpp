#ifndef MORT_COMMON_ATOMFUN_HPP
#define MORT_COMMON_ATOMFUN_HPP

#include <string>

namespace mort
{
    using std::string;

    class morf_t;

    static const int SINGLE_BOND = 1;
    static const int DOUBLE_BOND = 2;
    static const int TRIPLE_BOND = 3;
    static const int SP1 = 1;
    static const int SP2 = 2;
    static const int SP3 = 3;

    int get_atom_pie(const morf_t& a);

    int get_atom_hybrid(const morf_t& a);

    int get_atom_degree(const morf_t& a);

    int get_atom_valence(const morf_t& a);

    int get_atom_hydrogen(const morf_t& a);

    int count_single_bond(const morf_t& a);

    int count_double_bond(const morf_t& a);

    int count_triple_bond(const morf_t& a);

    bool is_sp1(const morf_t& a);

    bool is_sp2(const morf_t& a);

    bool is_sp3(const morf_t& a);

    bool is_amide_carbon( const morf_t& a );

    // get atom of a mo
    morf_t atom_1st(const morf_t& o);

    morf_t atom_2nd(const morf_t& o);

    morf_t atom_3rd(const morf_t& o);

    morf_t atom_4th(const morf_t& o);

    morf_t atom_oth(const morf_t& b, const morf_t& a);

    //
    string uniq_name( const morf_t& a );

    // a set of functions find, create atom for a residue
    //morf_t create_atom( morf_t& r, const string& name );

    //morf_t frcget_atom(morf_t& resd, const string& name);
    
    //morf_t get_atom(const morf_t& resd, const string& name);

    //bool get_atom(const morf_t& resd, const string& name, morf_t& atom);

 
} // namespace mort

#endif

