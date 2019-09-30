#include <sstream>
#include <boost/python.hpp>
#include <common.hpp>

using namespace mort;

using namespace boost::python;

void numvec_setitem(numvec& v, int id, double d)
{
    v[id] = d;
}

double numvec_getitem(numvec& v, int id)
{
    if( id < 0 || id >= v.size() )
    {
        throw std::out_of_range( "Error: out of range in numvec");
    }

    return v[id];
}

string numvec_tostr(const numvec& v)
{
    std::ostringstream os;
    os << v.size() << " ";
    os << "[";

    for( int i=0; i < v.size(); ++i)
    {
        os << v[i];
	os << ( i==v.size()-1 ? ']' : ',' );
    }
    return os.str();
}    

BOOST_PYTHON_MODULE(libpycommon)
{
    def( "hash",   mort::hash );
    def( "unhash", mort::unhash );
    
    class_<numvec>( "numvec", init< int >() )
        .def("size", &numvec::size )
        .def("__setitem__", &numvec_setitem )
	.def("__getitem__", &numvec_getitem )
	.def("__str__", &numvec_tostr )
	.def("__len__", &numvec::size )
    ;
}


