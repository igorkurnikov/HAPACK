#include <boost/python.hpp>
#include <object.hpp>
#include <capbox.hpp>
#include <capbox/addions.hpp>

using namespace mort;

using namespace boost::python;


struct solute_wrap : solute_i,wrapper<solute_i>
{
    int natom() const
    {
        return this->get_override("natom")();
    }

    vector<double>& getcord()
    {
        return this->get_override("getcord")();
    }
    
    vector<double>& getvdwr()
    {
        return this->get_override("getvdwr")();
    }

    void insert_resd( const molecule_t& svt, const numvec& sft, int idx )
    {
        this->get_override("insert_resd")( svt, sft, idx );
    }
  
    void finish()
    {
         this->get_override("finish")();
    }
};

struct ionee_wrap  : ionee_i, wrapper<ionee_i>
{
    int natom()
    {
        return this->get_override("natom")();
    }

    int nresd()
    {
        return this->get_override("nresd")();
    }

    int solute_nresd()
    {
        return this->get_override("solute_nresd")();
    }

    int solute_natom()
    {
        return this->get_override("solute_natom")();
    }

    vector<double>& getcord()
    {
        return this->get_override("getcord")();
    }

    vector<double>& getchrg()
    {
        return this->get_override("getchrg")();
    }

    vector<double>& getvdwr()
    {
        return this->get_override("getvdwr")();
    }

    vector<string>& gettype()
    {
        return this->get_override("gettype")();
    }

    void insert_resd( const molecule_t& ion, const numvec& pos, int idx)
    {
        this->get_override("insert_resd")( ion, pos, idx );
    }

    void remove_resd( int rid )
    {
        this->get_override("remove_resd")( rid );
    }

};

BOOST_PYTHON_MODULE(libpycapbox)
{
    class_<ionee_wrap, boost::noncopyable>("ionee_i")
        .def("natom", pure_virtual(&ionee_i::natom) )
        .def("nresd", pure_virtual(&ionee_i::nresd) )
        .def("solute_natom", pure_virtual(&ionee_i::solute_natom) )
        .def("solute_nresd", pure_virtual(&ionee_i::solute_nresd) )
        .def("getcord", pure_virtual(&ionee_i::getcord), return_internal_reference<1>() )
        .def("getchrg", pure_virtual(&ionee_i::getchrg), return_internal_reference<1>() )
        .def("getvdwr", pure_virtual(&ionee_i::getvdwr), return_internal_reference<1>() )
        .def("gettype", pure_virtual(&ionee_i::gettype), return_internal_reference<1>() )
        .def("insert_resd", pure_virtual(&ionee_i::insert_resd) )
        .def("remove_resd", pure_virtual(&ionee_i::remove_resd) ) 
       ;

    class_<solute_wrap, boost::noncopyable>("solute_i")
        .def("natom", pure_virtual(&solute_i::natom) )
        .def("getcord", pure_virtual(&solute_i::getcord), return_internal_reference<1>() )
        .def("getvdwr", pure_virtual(&solute_i::getvdwr), return_internal_reference<1>() )
        .def("insert_resd", pure_virtual(&solute_i::insert_resd) )
        .def("finish",  pure_virtual(&solute_i::finish) )
        ;
 
    def("_addions",      &addions);
    def("_addions_core", &addions_core);

    def("_solvatecap",   &solvatecap);
    def("_solvatebox",   &solvatebox);
    def("_solvateoct",   &solvateoct);
    def("_solvateshl",   &solvateshl);

    def("_solvatecap_core", &solvatecap_core);
    def("_solvatebox_core", &solvatebox_core);
    def("_solvateoct_core", &solvateoct_core);
    def("_solvateshl_core", &solvateshl_core);




}
