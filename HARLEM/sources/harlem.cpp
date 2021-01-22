#include <windows.h>
#include <iostream>
#include <vector>
#include <string>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

//! Start harlem using Harlem with python's wxApp
int start_harlem(int argc, char **argv)
{
	wchar_t *program = Py_DecodeLocale("HARLEM CONSOLE", NULL);

	Py_SetProgramName(program);  /* optional but recommended */
	Py_Initialize();

	std::vector<std::wstring> python_args;
	python_args.resize(argc);
	int i;
	wchar_t** argv_p = (wchar_t**)malloc(argc * sizeof(wchar_t*));
	for (i = 0; i < argc; i++)
	{
		std::string tmp(argv[i]);
		python_args[i].resize(tmp.length());
		std::copy(tmp.begin(), tmp.end(), python_args[i].begin());
		argv_p[i] = (wchar_t*)python_args[i].c_str();
	}

	PySys_SetArgv(argc, argv_p);

	PyRun_SimpleString("import molset.harlempy");
	PyRun_SimpleString("molset.harlempy.start_harlem()");

	//delete[] argv_p;

	if (Py_FinalizeEx() < 0) {
		exit(120);
	}
	PyMem_RawFree(program);
	return 0;
}

int main(int argc, char **argv)
{
	SetConsoleTitle(TEXT("HARLEM CONSOLE"));
	return start_harlem(argc, argv);
}
