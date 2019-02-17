#include <wx/wx.h>

int main(int argc, char **argv)
{
	SetConsoleTitle(TEXT("HARLEM CONSOLE"));
#ifdef _DEBUG
	LoadLibrary(TEXT("harlempy\\_halib_d.pyd"));
	LoadLibrary(TEXT("harlempy\\_molset_d.pyd"));
#else
	LoadLibrary(TEXT("harlempy\\_halib.pyd"));
	LoadLibrary(TEXT("harlempy\\_molset.pyd"));
#endif
	return wxEntry(argc, argv);
}
