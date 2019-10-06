#include <boost/python.hpp>
#include <mortgl.hpp>

using namespace mort;

using namespace boost::python;

struct graphic_wrapper : graphic_t, wrapper<graphic_t>
{
    string name() const
    {
        this->get_override("name")();
    }

    void paint()
    {
        this->get_override("paint")();
    }

    void update()
    {
        this->get_override("update")();
    }

    shared_ptr<graphic_t> clone(const molecule_t& m, const vector<string>& args) const
    {
        return this->get_override("clone")(m, args);
    }

};


BOOST_PYTHON_MODULE(libpymortgl)
{

    class_<drawing_t>( "drawing_t" )
        .def("add", &drawing_t::add)
	.def("has", &drawing_t::has)
	.def("get", &drawing_t::get)
	.def("remove", &drawing_t::remove)
	.def("mouse_press", &drawing_t::mouse_press)
	.def("mouse_move",  &drawing_t::mouse_move)
	.def("mouse_release", &drawing_t::mouse_release)
	.def("init", &drawing_t::init)
	.def("repaint", &drawing_t::repaint)
	.def("resize", &drawing_t::resize)
    ;

    class_<graphic_wrapper, boost::noncopyable>( "graphic_t" )
        .def("set", &graphic_t::set)
	.def("setall", &graphic_t::setall)
    ;


    objects::class_value_wrapper<
         shared_ptr<graphic_t>, 
	 objects::make_ptr_instance<graphic_t, objects::pointer_holder<shared_ptr<graphic_t>,graphic_t> >
	 >();


}


