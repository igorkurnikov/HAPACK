#include <sstream>
#include <boost/python.hpp>
#include <common.hpp>

using namespace mort;

using namespace boost::python;

template<typename vector_T, typename T >
void vector_setitem(vector_T& v, int id, T d)
{
    v[id] = d;
}

template<typename vector_T, typename T >
T vector_getitem(const vector_T& v, int id)
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
        .def("__setitem__", 	&vector_setitem<numvec, double> )
	.def("__getitem__", 	&vector_getitem<numvec, double> )
	.def("__str__", 	&numvec_tostr )
	.def("__len__", 	&numvec::size )
    ;

    class_< vector<double> >( "dbevec", init< int >() )
        .def("__setitem__", 	&vector_setitem<vector<double>, double> )
        .def("__getitem__", 	&vector_getitem<vector<double>, double> )
        .def("__len__", 	&vector<double>::size )
	.def("size", 		&vector<double>::size )
        .def("push_back", 	&vector<double>::push_back )
    ;

    class_< vector<int> >( "intvec", init< int >() )
        .def("__setitem__", 	&vector_setitem<vector<int>, int> )
        .def("__getitem__", 	&vector_getitem<vector<int>, int> )
        .def("__len__", 	&vector<int>::size )
        .def("size", 		&vector<int>::size )
        .def("push_back", 	&vector<double>::push_back )
      
    ;

    class_< vector<string> >( "strvec", init< int >() )
        .def("__setitem__", 	&vector_setitem<vector<string>, string> )
        .def("__getitem__", 	&vector_getitem<vector<string>, string> )
        .def("__len__", 	&vector<string>::size )
        .def("size", 		&vector<string>::size )
        .def("push_back", 	&vector<string>::push_back )
   ;


}


