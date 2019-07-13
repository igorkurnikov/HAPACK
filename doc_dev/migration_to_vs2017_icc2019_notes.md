# Migration to Visual Studio 2017 and Intel Compilers 2019

## Visual Studio

Visual Studio 2017 community edition is used
"https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=Community&rel=15#''

Installed packages
* Desktop development with C++ 
  ** C++/CLI support
  ** Visual C++ MFC for x86 and x64
* Python development 

## Intel Compilers

Intel® Parallel Studio XE Cluster Edition for Windows 2019 Update 1
It removed Intel Parallel Studio 2013, it can be later reinstalled.

## VCPKG



```
# install vcpkg
cd C:\PROG_SRC\
git clone https://github.com/nsimakov/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
# install libraries
.\vcpkg.exe install boost:x86-windows
.\vcpkg.exe install python2:x86-windows
.\vcpkg.exe install wxwidgets:x86-windows
.\vcpkg.exe install plplot[wxwidgets]:x86-windows
.\vcpkg.exe install mpir:x86-windows

```





## Code Changes and Problems During Migration

### AMBER11

AMBER11_IGOR git repo: https://gitlab.com/mkurnikovagroup/AMBER11_IGOR

C:\MYPROG\AMBER11_IGOR\AmberTools\src\gleap\mortsrc\common\hashcode.cpp
added #include <algorithm> for std:min

C:\MYPROG\AMBER11_IGOR\AmberTools\src\gleap\mortsrc\enefrc\nonbond-pbc.cpp
disable erfc declaration, it is in math.h

resd.cpp commented assert( is != NULL ), not sure does it have sense to compare istream& and NULL

### Intel Compilers and MKL

Intel Parallel Studio XE Cluster Edition for Windows 2019 Update 1, changed locations

C:\Program Files (x86)\Intel\Composer XE\compiler\lib\ia32

C:\Program Files (x86)\Intel\Composer XE\mkl\lib\ia32



Problem:
undefined reference to `for__rtc_uninit_use_src'
Resolution:
Setting Runtime checking to None (/check:none)
Custom and ensuring that Check Uninitialized Variables" is set to No
Ref: Remove the flag -check all,noarg_temp_created from FFLAGS and FCFLAGS
https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/601355



Problem:
On linking, complains about linking flag for SAFESEH:
2>mkl_intel_c_dll.lib(_dgesv.obj) : error LNK2026: module unsafe for SAFESEH image.
2>mkl_intel_c_dll.lib(_dsyev.obj) : error LNK2026: module unsafe for SAFESEH image.
2>mkl_intel_c_dll.lib(_dsygv.obj) : error LNK2026: module unsafe for SAFESEH image.
2>mkl_intel_c_dll.lib(_dgetrf.obj) : error LNK2026: module unsafe for SAFESEH image.
2>mkl_intel_c_dll.lib(_dgetri.obj) : error LNK2026: module unsafe for SAFESEH image.
Resolution:
SAFESEH:NO

Problem (in Debug):

1>gmathc.obj : error LNK2001: unresolved external symbol __imp____argc
1>gmathc.obj : error LNK2001: unresolved external symbol __imp____argv

Resolution:

Alternative implementation of gmathc.iargc_ and getarg_ taken from jumna to tackle 'unresolved external symbol __imp____argc' 



### Boost

Boost: version 1.68.0
boost_filesystem-vc140-mt.lib need to be specified directly


### GMP

mpir used instead

https://github.com/Microsoft/vcpkg/issues/2628

## wxWidgets

debug version of wxWidgets do asserts on absence of  vertical alignment for wxBoxSizer( wxVERTICAL ) and similar for horizontal, e.g.:

load_dlg.sizer_main_v->Add( calc_bonds_chk, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

should be

load_dlg.sizer_main_v->Add( calc_bonds_chk, 0, wxALL, 5 );

Resolution:

Disable checking in wxWidget library (patch is in our vcpkg repo)







# x64 bit

## VCPKG with boost and mpir for 64-bit
```cmd
# install vcpkg
cd C:\PROG_SRC\
# if vcpgk is not installed overwise use one installed with x86
#git clone https://github.com/nsimakov/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
# install libraries
.\vcpkg.exe install boost:x64-windows
.\vcpkg.exe install mpir:x64-windows

```










