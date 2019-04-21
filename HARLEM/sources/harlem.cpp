//#include <windows.h>
#include <wx/wx.h>

typedef int(__stdcall *f_start_harlemappwx)(int argc, char **argv);

int main(int argc, char **argv)
{
	SetConsoleTitle(TEXT("HARLEM CONSOLE"));
#ifdef _DEBUG
	HINSTANCE hLLPNPSDLL = LoadLibrary(TEXT("harlempy\\_llpnps_d.pyd"));
	HINSTANCE hMolSetDLL = LoadLibrary(TEXT("harlempy\\_molset_d.pyd"));
#else
	HINSTANCE hLLPNPSDLL = LoadLibrary(TEXT("harlempy\\_llpnps.pyd"));
	HINSTANCE hMolSetDLL = LoadLibrary(TEXT("harlempy\\_molset.pyd"));
#endif
	if (!hLLPNPSDLL) {
		std::cout << "Cannot load harlempy\\_llpnps.pyd" << std::endl;
		std::cout << "Last Error: " << GetLastError() << std::endl;
		Sleep(3);
		return EXIT_FAILURE;
	}
	if (!hMolSetDLL) {
		std::cout << "Cannot load harlempy\\_molset.pyd" << std::endl;
		std::cout << "Last Error: " << GetLastError() << std::endl;
		Sleep(3);
		return EXIT_FAILURE;
	}

	int status = wxEntry(argc, argv);

	FreeLibrary(hMolSetDLL);

	return status;
}
