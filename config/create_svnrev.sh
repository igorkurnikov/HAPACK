#!/usr/bin/sh
 echo '#include "hasvnrev.h"'> hasvnrev.cpp
 echo -n 'const char* HaSVNRevision(){const char* SVN_Version = "' >> hasvnrev.cpp
 echo -n `svn info ./|grep "Revision"` >> hasvnrev.cpp
 echo '"; return SVN_Version; }'   >> hasvnrev.cpp
 echo -n 'const char* HaSVNDate(){const char* SVN_Date = "' >> hasvnrev.cpp
 echo -n `svn info ./|grep "Last Changed Date"` >> hasvnrev.cpp
 echo '"; return SVN_Date; }'   >> hasvnrev.cpp
