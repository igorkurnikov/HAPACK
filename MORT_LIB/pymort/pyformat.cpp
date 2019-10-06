
#include <string>
#include <boost/python.hpp>
#include <object.hpp>
#include <format.hpp>


using namespace mort;

using namespace boost::python;

BOOST_PYTHON_MODULE(libpyformat)
{

    void (*load_mol_ptr)(const string&, molecule_t&, const hashid_t&) = &load_mol;
    void (*save_mol_ptr)(const string&, const molecule_t&, const hashid_t&) = &save_mol;
    void (*load_mdb_ptr)(const string&, database_t&, const hashid_t&) = &load_mdb;
    void (*save_mdb_ptr)(const string&, const database_t&, const hashid_t&) = &save_mdb;
    def("_load_mol", load_mol_ptr);
    def("_save_mol", save_mol_ptr);
    def("_load_mdb", load_mdb_ptr);
    def("_save_mdb", save_mdb_ptr);
    def("_load_frc",  &load_frc);

    def("save_top",   &save_top);
}
