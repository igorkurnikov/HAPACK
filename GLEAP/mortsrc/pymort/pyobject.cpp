#include <stdexcept>
#include <boost/python.hpp>
#include <common.hpp>
#include <object.hpp>
#include <pdbent.hpp>
using namespace mort;

using std::runtime_error;

using namespace boost::python;

template< typename T >
void set_parm(T* p, const string& parmname, const object& parmvalue)
{
    if( parmname.length()>=2 && parmname.substr(0,2) == "__" )
    {
        throw runtime_error( "Error: can not set built in attribute" );
    }

    extract<int> igetter(parmvalue);
    if( igetter.check() )
    {
        int ivalue = igetter();
	p->set_i(parmname, ivalue);
	return;
    }

    extract<double> dgetter(parmvalue);
    if( dgetter.check() )
    {
        double dvalue = dgetter();
	p->set_d(parmname, dvalue);
	return;
    }

    extract<string> sgetter(parmvalue);
    if( sgetter.check() )
    {
        string svalue = sgetter();
	p->set_s(parmname, svalue);
	return;
    }

    extract<numvec> vgetter(parmvalue);
    if( vgetter.check() )
    {
        numvec vvalue = vgetter();
	p->set_v(parmname, vvalue);
	return;
    }

/*
    any avalue = pyobj_to_any(o);
    p->set_a(parmname, avalue);
 */

}

template<typename T>
object get_parm(const T* p, const string& parmname)
{
    int ivalue;
    if( p->get_i(parmname, ivalue) )
    {
        return object(ivalue);
    }

    double dvalue;
    if( p->get_d(parmname, dvalue) )
    {
        return object(dvalue);
    }

    string svalue;
    if( p->get_s(parmname, svalue) )
    {
        return object(svalue);
    }

    numvec vvalue;
    if( p->get_v(parmname, vvalue) )
    {
        return object(vvalue);
    }

/*
    any avalue;
    if( p->get_a(parmname, avalue) )
    {
        return any_to_pyobj(avalue);
    }
 */

    throw runtime_error( "Error: entity do not have parm: " + parmname );
}

morf_t range_getitem(const mobj_range& r, int id)
{
    if( id<0 || id>=r.size() )
    {
        PyErr_SetString(PyExc_IndexError, "out of range");
        throw_error_already_set();
    }

    return r[id];
}

morf_t atomvec_getitem(const atomvec_t& r, int id)
{
    if( id<0 || id>=r.size() )
    {
        PyErr_SetString(PyExc_IndexError, "out of range");
        throw_error_already_set();
    }

    return r[id];
}


python::list database_keys(const database_t& db)
{
   python::list keys;
   database_t::const_iterator i = db.begin();
   for( ; i != db.end(); ++i )
   {
       keys.append( i->first );
   }

   return keys;
}


BOOST_PYTHON_MODULE(libpyobject)
{
    class_<root_t>( "root_t" )
        .def("__setattr__", &set_parm<root_t>)
        .def("__getattr__", &get_parm<root_t>)
    ;

    objects::class_value_wrapper<
         shared_ptr<root_t>, 
	 objects::make_ptr_instance<root_t, objects::pointer_holder<shared_ptr<root_t>,root_t> >
	 >();


    class_<molecule_t, bases<root_t> >( "molecule_t" )
        .add_property("natom", &molecule_t::natom)
	.add_property("nbond", &molecule_t::nbond)
	.add_property("nangl", &molecule_t::nangl)
	.add_property("ntor2", &molecule_t::ntor2)
	.add_property("nptor", &molecule_t::nptor)
	.add_property("nresd", &molecule_t::nresd)
	.add_property("atoms", &molecule_t::atoms)
	.add_property("bonds", &molecule_t::bonds)
	.add_property("angls", &molecule_t::angls)
	.add_property("dihes", &molecule_t::dihes)
	.add_property("tor2s", &molecule_t::tor2s)
	.add_property("imprs", &molecule_t::imprs)
	.add_property("ptors", &molecule_t::ptors)
	.add_property("resds", &molecule_t::resds)
    ;
    
    objects::class_value_wrapper<
         shared_ptr<molecule_t>, 
	 objects::make_ptr_instance<molecule_t, objects::pointer_holder<shared_ptr<molecule_t>,molecule_t> >
	 >();


    class_<mobj_range>("morange_t", init<const mobjiter_t&, const mobjiter_t&>() )
         .def("__getitem__", &range_getitem)
	 .def("__len__", &mobj_range::size)
	 .def("size", &mobj_range::size)
    ;

    class_<morf_t>( "morf_t", init< const molecule_t&, const hashid_t&, int >() )
        .def("__setattr__", &set_parm<morf_t>)
	.def("__getattr__", &get_parm<morf_t>)
	.def("connect", &morf_t::connect)
	.def("disconnect", &morf_t::disconnect)
	.def("is_connected_to", &morf_t::is_connected_to)
	.add_property("atoms", &morf_t::atoms)
	.add_property("bonds", &morf_t::bonds)
	.add_property("angls", &morf_t::angls)
	.add_property("resds", &morf_t::resds)
	.add_property("absid", &morf_t::absid)
	.add_property("relid", &morf_t::relid)
    ;

    class_<database_t, bases<root_t> >( "database_t" )
	.def("has", &database_t::has)
        .def("set", &database_t::set)
	.def("get", &database_t::get) //  , return_internal_reference<1>() )
	.def("has_key", &database_t::has)
        .def("__setitem__", &database_t::set)
        .def("__getitem__", &database_t::get)
        .def("get_mol", &database_t::get_mol)
	.def("get_mdb", &database_t::get_mdb)
	.def("keys", &database_keys)
    ;

    class_<atomvec_t, bases<root_t> >( "atomvec_t" )
        .def("__getitem__", &atomvec_getitem)
	.def("__len__", &atomvec_t::size)
	.def("size", &atomvec_t::size)
    ;


    class_<namemap_t, bases<root_t> >( "namemap_t" )
        .def( "add_resd_map", &namemap_t::add_resd_map )
        .def( "add_atom_map", &namemap_t::add_atom_map )
    ;

}


