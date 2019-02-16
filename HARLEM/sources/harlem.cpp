#include <wx/wx.h>

int main(int argc, char **argv)
{
	SetConsoleTitle(TEXT("HARLEM CONSOLE"));
#ifdef _DEBUG
	LoadLibrary(TEXT("_molset_d.pyd"));
#else
	LoadLibrary(TEXT("_molset.pyd"));
#endif
	return wxEntry(argc, argv);
}
