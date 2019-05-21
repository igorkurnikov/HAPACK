#include <windows.h>
#include <iostream>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

typedef int(__stdcall *f_start_harlemappwx)(int argc, char **argv);

//! Start harlem using HarlemAppWX
int start_harlemappwx(int argc, char **argv)
{
#ifdef _DEBUG
	//HINSTANCE hLLPNPSDLL = LoadLibrary(TEXT("harlempy\\_llpnps_d.pyd"));
	HINSTANCE hMolSetDLL = LoadLibrary(TEXT("harlempy\\_molset_d.pyd"));
#else
	//HINSTANCE hLLPNPSDLL = LoadLibrary(TEXT("harlempy\\_llpnps.pyd"));
	HINSTANCE hMolSetDLL = LoadLibrary(TEXT("harlempy\\_molset.pyd"));
#endif
	/*if (!hLLPNPSDLL) {
		std::cout << "Cannot load harlempy\\_llpnps.pyd" << std::endl;
		std::cout << "Last Error: " << GetLastError() << std::endl;
		Sleep(3000000);
		return EXIT_FAILURE;
	}*/
	if (!hMolSetDLL) {
		std::cout << "Cannot load harlempy\\_molset.pyd" << std::endl;
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

//! Start harlem using Harlem with python's wxApp
int start_harlem(int argc, char **argv)
{
	wchar_t *program = Py_DecodeLocale("HARLEM CONSOLE", NULL);

	Py_SetProgramName(program);  /* optional but recommended */
	Py_Initialize();

	PyRun_SimpleString(
		"from harlem.start_harlem import start_harlem\n"
		"start_harlem()\n");

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
