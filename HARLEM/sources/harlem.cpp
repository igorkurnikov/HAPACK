#include <windows.h>
#include <iostream>
#include <vector>
#include <string>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#if defined(_MSC_VER)
typedef int(__stdcall *f_start_harlemappwx)(int argc, char **argv);

//! Start harlem using HarlemAppWX
int start_harlemappwx(int argc, char **argv)
{
#ifdef _DEBUG
	//HINSTANCE hLLPNPSDLL = LoadLibrary(TEXT("molsetll\\_llpnps_d.pyd"));
	HINSTANCE hMolSetDLL = LoadLibrary(TEXT("molsetll\\_molset_d.pyd"));
#else
	//HINSTANCE hLLPNPSDLL = LoadLibrary(TEXT("molsetll\\_llpnps.pyd"));
	HINSTANCE hMolSetDLL = LoadLibrary(TEXT("molset\\_molset.pyd"));
#endif
	/*if (!hLLPNPSDLL) {
		std::cout << "Cannot load molset\\_llpnps.pyd" << std::endl;
		std::cout << "Last Error: " << GetLastError() << std::endl;
		Sleep(3000000);
		return EXIT_FAILURE;
	}*/
	if (!hMolSetDLL) {
		std::cout << "Cannot load molset\\_molset.pyd" << std::endl;
		std::cout << "Last Error: " << GetLastError() << std::endl;
		Sleep(3000000);
		return EXIT_FAILURE;
	}

#if defined(_M_IX86)
	f_start_harlemappwx start_harlemappwx = (f_start_harlemappwx)GetProcAddress(hMolSetDLL, "_start_harlemappwx@8");
#else
	f_start_harlemappwx start_harlemappwx = (f_start_harlemappwx)GetProcAddress(hMolSetDLL, "start_harlemappwx");
#endif
	if (!start_harlemappwx) {
		std::cout << "could not locate the function" << std::endl;
		Sleep(3000000);
		return EXIT_FAILURE;
	}

	//int status = wxEntry(argc, argv);
	int status = start_harlemappwx(argc, argv);

	FreeLibrary(hMolSetDLL);
	return status;
}
#endif

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

	PyRun_SimpleString(
		"from harlempy import *\n"
		"from harlempy.start_harlem import start_harlem\n"
		"start_harlem()\n");

	delete[] argv_p;

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
