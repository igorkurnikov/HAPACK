#include <boost/python.hpp>
#include <common.hpp>
#include <object.hpp>
#include <format.hpp>
#include <guilib.hpp>
#include <capbox.hpp>
#include <atmask.hpp>

using namespace mort;

using namespace boost::python;

template< typename T >
void set_parm(T* p, const string& parmname, const object& parmvalue)
{
    if( parmname.substr(0,2) == "__" )
    {
        throw logic_error( "Error: can not set built in attribute" );
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

    throw logic_error( "Error: entity do not have parm: " + parmname );
}

moref_t range_getitem(const range_t& r, int id)
{
    if( id<0 || id>=r.size() )
    {
        PyErr_SetString(PyExc_IndexError, "out of range");
        throw_error_already_set();
    }

    return r[id];
}

moref_t atomvec_getitem(const atomvec_t& r, int id)
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
    class_<entity_t>( "entity_t" )
        .def("__setattr__", &set_parm<entity_t>)
        .def("__getattr__", &get_parm<entity_t>)
    ;

    objects::class_value_wrapper<
         shared_ptr<entity_t>, 
	 objects::make_ptr_instance<entity_t, objects::pointer_holder<shared_ptr<entity_t>,entity_t> >
	 >();


    class_<molecule_t, bases<entity_t> >( "molecule_t" )
	.def("create_atom", &molecule_t::create_atom)
	.def("create_bond", &molecule_t::create_bond)
	.def("create_angl", &molecule_t::create_angl)
	.def("create_resd", &molecule_t::create_resd)
        .add_property("natom", &molecule_t::natom)
	.add_property("nbond", &molecule_t::nbond)
	.add_property("nangl", &molecule_t::nangl)
	.add_property("nresd", &molecule_t::nresd)
	.add_property("atoms", &molecule_t::atoms)
	.add_property("bonds", &molecule_t::bonds)
	.add_property("angls", &molecule_t::angls)
	.add_property("resds", &molecule_t::resds)
    ;
    
    objects::class_value_wrapper<
         shared_ptr<molecule_t>, 
	 objects::make_ptr_instance<molecule_t, objects::pointer_holder<shared_ptr<molecule_t>,molecule_t> >
	 >();


    class_<range_t>("range_t", init<const moiter_t&, const moiter_t&>() )
         .def("__getitem__", &range_getitem)
	 .def("__len__", &range_t::size)
	 .def("size", &range_t::size)
    ;

    class_<moref_t>( "moref_t", init< const molecule_t&, const hashid_t&, int >() )
        .def("__setattr__", &set_parm<moref_t>)
	.def("__getattr__", &get_parm<moref_t>)
	.def("connect", &moref_t::connect)
	.def("disconnect", &moref_t::disconnect)
	.def("is_connected_to", &moref_t::is_connected_to)
	.add_property("related_atoms", &moref_t::related_atoms)
	.add_property("related_bonds", &moref_t::related_bonds)
	.add_property("related_angls", &moref_t::related_angls)
	.add_property("related_resds", &moref_t::related_resds)
    ;

    class_<database_t, bases<entity_t> >( "database_t" )
        .def("set", &database_t::set)
	.def("has", &database_t::has)
	.def("get", &database_t::get) //  , return_internal_reference<1>() )
	.def("get_mol", &database_t::get_mol)
	.def("get_mdb", &database_t::get_mdb)
	.def("keys", &database_keys)
    ;

    void (*load_mol_fn)(const string&, molecule_t&, int) = &load_mol;
    void (*save_mol_fn)(const string&, const molecule_t&, int) = &save_mol;
    void (*load_mdb_fn)(const string&, database_t&, int) = &load_mdb;
    void (*save_mdb_fn)(const string&, const database_t&, int) = &save_mdb;
    atomvec_t (*ambmsk_mask_atom_fn)(const molecule_t&, const string& ) = &mask_atom;
    atomvec_t (*smarts_mask_atom_fn)(const molecule_t&, const string& ) = &smarts_mask_atom;

    def("create_bond", &create_bond);

    def("_ambmsk_mask_atom", ambmsk_mask_atom_fn);
    def("_smarts_mask_atom", smarts_mask_atom_fn);
    def("_load_mol", load_mol_fn);
    def("_save_mol", save_mol_fn);
    def("_load_mdb", load_mdb_fn);
    def("_save_mdb", save_mdb_fn);
  
    double (*molcharge)(const molecule_t&) = &charge;
    double (*rescharge)(const moref_t& ) = &charge;

    def("_molcharge", molcharge);
    def("_rescharge", rescharge);

    def("add_ion",     &add_ion);
    def("solvate_cap", &solvate_cap);
    def("solvate_box", &solvate_box);
    def("solvate_oct", &solvate_oct);
    def("solvate_shell", &solvate_shell);
/*
    def("_load_amber_ffp",  &load_amber_ffp);
    def("_load_amoeba_ffp", &load_amoeba_ffp);
    def("assign_ffparm", &assign_ffparm);
 */

    class_<atomvec_t, bases<entity_t> >( "atomvec_t" )
        .def("__getitem__", &atomvec_getitem)
	.def("__len__", &atomvec_t::size)
	.def("size", &atomvec_t::size)
    ;

/*
    def("eval_bond", &eval_bond);
    def("eval_angl", &eval_angl);
    def("eval_tors", &eval_tors);
    def("eval_oops", &eval_oops);
    def("eval_direct", &eval_direct);
    def("eval_regewald", &eval_regewald);
 */
}


