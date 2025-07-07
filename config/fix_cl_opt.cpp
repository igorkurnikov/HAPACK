#include <string>
#include <vector>
#include <process.h>
#include <errno.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	std::string cmd_line;
	std::string exe_name;
	std::string tmp;
	std::string sa;

	std::vector<std::string>arg_vec;
	
	int i,j;
	char dquote = '\"';
	
	enum FORTRAN_TYPE { COMPAQ_FORTRAN_V6 = 1, INTEL_FORTRAN_IA32_V11_072 = 2 }; 
	enum CPP_VERSION { VISUAL_CPP_8 = 1 }; 
	
	int fortran_ver = INTEL_FORTRAN_IA32_V11_072;
//	int fortran_ver = COMPAQ_FORTRAN_V6;
	int cpp_ver     = VISUAL_CPP_8;
	
	std::string vs_install_dir;
	std::string vc_install_dir;
	std::string framework_dir;
	std::string framework_version;
	std::string framework_sdk_dir;
	std::string fortran_install_dir;
	
	std::string compaq_fortran_v6_install_dir;
		
#if defined(_MSC_VER)
	if( cpp_ver == VISUAL_CPP_8)
	{
		vs_install_dir = "C:\\Program Files\\Microsoft Visual Studio 8";
		vc_install_dir = "C:\\Program Files\\Microsoft Visual Studio 8\\VC";
		framework_dir = "C:\\WINDOWS\\Microsoft.NET\\Framework";
		framework_version = "v2.0.50727";
		framework_sdk_dir = "C:\\Program Files\\Microsoft Visual Studio 8\\SDK\\v2.0";
	}
	compaq_fortran_v6_install_dir = "C:\\Program Files\\Microsoft Visual Studio\\DF98";
	if( fortran_ver == COMPAQ_FORTRAN_V6 )
	{
	     fortran_install_dir= compaq_fortran_v6_install_dir;
	}
	else if( fortran_ver == INTEL_FORTRAN_IA32_V11_072 )
	{
             fortran_install_dir= "C:\\Program Files\\Intel\\Compiler\\11.0\\072\\fortran";
	}

#else
	if( cpp_ver == VISUAL_CPP_8)
	{
		 vs_install_dir = "/cygdrive/C/Program Files/Microsoft Visual Studio 8";
	         vc_install_dir = "/cygdrive/C/Program Files/Microsoft Visual Studio 8/VC";
	         framework_dir = "/cygdrive/C/WINDOWS/Microsoft.NET/Framework";
	         framework_version = "v2.0.50727";
	         framework_sdk_dir = "/cygdrive/c/Program Files/Microsoft Visual Studio 8/SDK/v2.0";
	}
	compaq_fortran_v6_install_dir = "/cygdrive/c/Program Files/Microsoft Visual Studio/DF98";
	if( fortran_ver == COMPAQ_FORTRAN_V6 )
	{
	     fortran_install_dir= compaq_fortran_v6_install_dir;
	}
	else if( fortran_ver == INTEL_FORTRAN_IA32_V11_072 )
	{
             fortran_install_dir= "/cygdrive/c/Program Files/Intel/Compiler/11.0/072/fortran";
	}
#endif
	
	std::string str;
	std::string mpi_libs;
	std::string mpi_libs_path;
	
	mpi_libs = "mpi.lib fmpich2h.lib";
	mpi_libs_path = "C:\\Program Files\\MPICH2\\LIB";

	int len_sa;
	char* ptr;
	int nl;
	
	int has_space = 0;
	
	int dash_c   = 0;
	std::string arg_o;  // value of -o option
		
	int is_link= 0;
	int is_fortran= 0;
	int is_debug= 0;
	
	for(i=1; i < argc; i++)
	{
		sa = argv[i];
		len_sa = sa.length();
		if( sa == "-vs_link"  )
		{
			is_link = 1;
			continue;
		}
		if( sa == "-fortran"  )
		{
			is_fortran = 1;
			continue;
		}
		if( sa == "-fortran_dir"  )
		{
			i++;
			fortran_install_dir = argv[i];
			sa = "";
			continue;
		}

		has_space = 0;
		for(j=0; j < sa.length(); j++)
		{
			if(sa[j] == 0x20) has_space = 1;
		}
		if(has_space) 
		{ 
			sa = dquote + sa;
			sa += dquote;
		}
		if(sa == "-all-static") sa = "";
		if(sa == "--version") sa = "";
		if(sa == "-Wno-non-template-friend") sa = "";
		if(sa == "-Wno-deprecated") sa = "";
		if(sa == "-ffixed-line-length-132") sa = "/extend_source:132";
		if(sa == "-O3") sa = "/Od";
		if(sa == "-O0") sa = "/Od";
		if(sa == "-O1") sa = "/Od";
		if(sa == "-O2") sa = "/Od";
		if(sa == "-Os") sa = "/Od"; // optimize for size
		
		if(sa.substr(0,2) == "-L")
		{
			if(is_link)
			{
				sa = "/libpath:" + sa.substr(2);
			}
			else
			{
				sa = "/link /libpath:" + sa.substr(2);
			}
		}
		
		if(sa == "-g")
		{
			is_debug = 1;
			continue;
		}
		
		if(sa == "-o") 
		{ 
			sa = "";
			i++;
			if( i >= argc) break;
			arg_o = argv[i]; 
		}
		if(sa == "-c") { dash_c = 1;  sa = "/c"; }
		if(sa == "/c") { dash_c = 1;  sa = "/c"; }
		
		if(sa == "/MD" && is_link) sa ="";
		if(sa == "/vd0" && is_link) sa ="";
		if(sa == "/EHsc" && is_link) sa ="";
		if((sa == "/Ox" || sa == "/Od" || sa == "/O1" || sa == "/O2") && is_link) sa ="";
		if(sa == "/Ox") sa ="/Od";
		
		if(sa.length() > 0) 
		{ 
			arg_vec.push_back(sa); 
			cmd_line += " "; cmd_line += sa;
		}
	}
	
	if(  !is_fortran && !is_link)
	{
	     if(is_debug)
	     {
		  cmd_line = " /nologo /MD /vd0 /EHs /Od /Z7 " + cmd_line; 
	     }
	     else
	     {
		  cmd_line = " /nologo /MD /vd0 /EHs /Ox " + cmd_line;
	     }
	}
	else if( is_fortran && !is_link)
	{
	     if(is_debug)
	     {
		  cmd_line = " /nologo /Z7 /assume:underscore /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /threads /Ox " + cmd_line; 
	     }
	     else
	     {
		  cmd_line = " /nologo /assume:underscore /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /threads " + cmd_line;
	     }
	}
	else if( is_link )
	{
	     if(is_debug)
	     {
		  cmd_line = " /debug " + cmd_line;
	     }
	     else
	     {
		  
	     }
	     str = "\""; str += "/libpath:";  
	     cmd_line = str + mpi_libs_path + "\" " + mpi_libs + cmd_line;
	}
	
		
	if( arg_o.length() > 0 )
	{
		if( dash_c )
		{
			tmp="/Fo";
		}
		else
		{
			tmp="/Fe";
			if(is_link) tmp="/out:";
		}
		cmd_line +=" ";
		cmd_line += tmp + arg_o;
		arg_o = "";
	}	

	std::string path_str = "";
	std::string include_str = "";
	std::string lib_str = "";
	std::string libpath_str = "";
	
	path_str += vs_install_dir + "\\Common7\\IDE;";
	path_str += vc_install_dir + "\\BIN;";
	path_str += vs_install_dir + "\\Common7\\Tools\\bin:";
	path_str += vc_install_dir + "\\PlatformSDK\\bin:";
	path_str += framework_sdk_dir + "\\bin:";
	path_str += framework_dir + "\\" + framework_version + ":";
	path_str += vc_install_dir + "\\VCPackages:";
//	path_str += getenv("PATH");
	
	include_str += vc_install_dir + "\\ATLMFC\\INCLUDE;";
	include_str += vc_install_dir + "\\INCLUDE;";
	include_str += vc_install_dir + "\\PlatformSDK\\include;";
	include_str += framework_sdk_dir + "\\INCLUDE;";
	
	lib_str  += vc_install_dir + "\\ATLMFC\\LIB;";
	lib_str  += vc_install_dir + "\\LIB;"; 
	lib_str  += vc_install_dir + "\\PlatformSDK\\lib;"; 
	lib_str  += framework_sdk_dir + "\\lib;"; 
	
	libpath_str += framework_dir + "\\" + framework_version +  ";";
	libpath_str += vc_install_dir + "\\ATLMFC\\LIB";

	if( fortran_ver == COMPAQ_FORTRAN_V6 )
	{
		path_str    = fortran_install_dir + "\\BIN;" + path_str;
		include_str = fortran_install_dir + "\\INCLUDE;" + include_str;
		lib_str       = fortran_install_dir + "\\LIB;" + lib_str;
	}
	else if ( fortran_ver == INTEL_FORTRAN_IA32_V11_072 )
	{
		path_str    = fortran_install_dir + "\\BIN\\IA32;" + path_str;
		include_str = fortran_install_dir + "\\INCLUDE;" + fortran_install_dir + "\\INCLUDE\\IA32;"  + include_str;
		lib_str       = fortran_install_dir + "\\LIB\\IA32;" + lib_str;
//		lib_str       = compaq_fortran_v6_install_dir + "\\LIB;" + lib_str;
		libpath_str = fortran_install_dir + "\\LIB\\IA32;" + libpath_str;
//		libpath_str = compaq_fortran_v6_install_dir + "\\LIB;" + libpath_str;
	}
	
	std::string intel_license_dir;
	std::string intel_license_str;
#if defined (_MSC_VER) 
	intel_license_dir = "C:\\Program Files\\Common Files\\Intel\\Licenses";
#else
	intel_license_dir = "/cygdrive/C/Program Files/Common Files/Intel/Licenses";
#endif
	
#if defined(_MSC_VER)
	_putenv_s("PATH",path_str.c_str());
	_putenv_s("INCLUDE",include_str.c_str());
	_putenv_s("LIB",lib_str.c_str());
	_putenv_s("LIBPATH",libpath_str.c_str());
	_putenv_s("INTEL_LICENSE_FILE",intel_license_dir.c_str());
#else
	path_str= "PATH=" + path_str;             putenv(path_str.c_str());
	include_str= "INCLUDE=" + include_str;  putenv(include_str.c_str());
	lib_str= "LIB=" + lib_str;                      putenv(lib_str.c_str());
	libpath_str= "LIBPATH=" + libpath_str;   putenv(libpath_str.c_str());
	intel_license_str = "INTEL_LICENSE_FILE=" + intel_license_dir; putenv(intel_license_str.c_str());
#endif
	
	int ires = 0;
	exe_name = "cl";

	if(is_link) exe_name = "link";
	if(is_fortran) 
	{
		if( fortran_ver == COMPAQ_FORTRAN_V6 )
		{
			exe_name = "df";
		}
		else if( fortran_ver == INTEL_FORTRAN_IA32_V11_072 )
		{
			exe_name = "ifort";
		}
	}

	
	std::string full_cmd_line = exe_name + "  ";
	full_cmd_line += cmd_line;
	
	printf("\nCMD LINE: %s\n",full_cmd_line.c_str());
	
#if defined(_MSC_VER)
	ires = _spawnlp(_P_WAIT,exe_name.c_str(),exe_name.c_str(),cmd_line.c_str(),NULL);
#else
	ires = spawnlp(P_WAIT,exe_name.c_str(),exe_name.c_str(),cmd_line.c_str(),NULL);
#endif
//	ires = spawnv(P_WAIT,exe_name.c_str(),exec_args);
	if(errno != 0)
	{
		printf("%s\n",sys_errlist[errno]);
	}
	printf("\n exec return code: %d \n",ires);
	return ires;
}