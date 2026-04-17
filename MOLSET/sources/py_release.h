#ifndef HARLEM_PY_RELEASE_H
#define HARLEM_PY_RELEASE_H

#ifdef _WIN32
#pragma push_macro("_DEBUG")
#undef _DEBUG
#include <Python.h>
#pragma pop_macro("_DEBUG")
#else
#include <Python.h>
#endif

#endif
