#include <boost/python.hpp>
#include <common.hpp>
#include <object.hpp>
#include <guilib.hpp>
#include <pdbent.hpp>

using namespace mort;

using namespace boost::python;


BOOST_PYTHON_MODULE(libpyobjfun)
{

/*
    def("create_bond", &create_bond);
    def("create_angl", &create_angl);
    def("create_tors", &create_tors);
    def("create_tor2", &create_tor2);
    def("create_ptor", &create_ptor);
*/


    double (*molcharge)(const molecule_t&) = &charge;
    double (*objcharge)(const morf_t& ) = &charge;

    def("_molcharge", molcharge);
    def("_rescharge", objcharge);

    def( "_mortenv", &mortenv, return_value_policy<reference_existing_object>() );
    def( "get_vdwr", &get_vdwr );

    def( "atom_1st", &atom_1st );
    def( "atom_2nd", &atom_2nd );
    def( "atom_3rd", &atom_3rd );
    def( "atom_4th", &atom_4th );

    def( "find_file", &find_file );
    def( "mdlize_mdb", &mdlize_mdb );
    def( "mdlize_seq", &mdlize_seq );

}
