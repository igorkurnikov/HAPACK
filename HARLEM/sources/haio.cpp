/*! haio.cpp 
  
   Classes and functions dealing with input and output 
   in HARLEM
 
   \author Igor Kurnikov
   \date 1998-2002
  
*/

#define HAIO_CPP

#include <mpi.h>

#include <wx/string.h>
#include <wx/log.h>

#include "haconst.h"
#include "stdarg.h"
#include "hastring.h"
#include "haio.h"

#if !defined(_MSC_VER)
#include <unistd.h>
#else 
#include <io.h>
#include <fcntl.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winbase.h>
#include <errno.h>
#endif

#include "hampi.h"
#include "harlemapp.h"


// extern "C" DllExport 
int PrintLog(const char* str, ... )
{
	va_list arg_list;
	va_start( arg_list, str );     /* Initialize variable arguments. */

	if( pApp != NULL && pApp->file_log != NULL)
	{
		vfprintf(pApp->file_log,str, arg_list);
		fflush(pApp->file_log);
	}
	else
	{
#if defined(HA_NOGUI) 
		vprintf(str, arg_list);
#else
		if( pApp == NULL )
		{
			vprintf(str, arg_list);
		}
		else if( pApp->mpi_driver != NULL && pApp->mpi_driver->myrank == 0 )
		{
	//		wxVLogGeneric(wxLOG_Message,str,arg_list);
	//		wxVLogMessage(str,arg_list);
			vprintf(str, arg_list);
		}
		else
		{
			vprintf(str, arg_list);
		}					
#endif
	}

//	wxLog::OnLog(1,wxString::FormatV(str, arg_list),0);
	va_end(arg_list);              /* Reset variable arguments.      */
	return TRUE;
}

 
void write_log_(const char* str, int n)
{
	int imax = n;
	for(; imax >=0; imax--)
	{
		if( !isspace(str[imax]) ) break;
	}

	if( imax > 0 )
	{
		int i;
		for( i = 0; i < imax; i++)
			PrintLog("%c",str[i]);
		PrintLog("\n");
	}
}

int ErrorMessage(const char* str)
{
	PrintMessage(str);
	return 1;
}

static FILE* con_out_fp=NULL;
static FILE* curr_stdout_fp=NULL;

void AddProcessNumberToFileName(char * out, const char *in, const char *pref, int pnum, int totnum)
{
	char in_main[512];
	char *in_ext, *in_det;
	strcpy(in_main, in);
	in_det = strrchr(in_main, '.');
	if (in_det != NULL)
	{
		*in_det = '\0';
		in_ext = in_det + 1;
		if (totnum <= 10)
			sprintf(out, "%s%s%.1d.%s\0", in_main, pref, pnum, in_ext);
		else if (totnum <= 100)
			sprintf(out, "%s%s%.2d.%s\0", in_main, pref, pnum, in_ext);
		else
			sprintf(out, "%s%s%.3d.%s\0", in_main, pref, pnum, in_ext);
	}
	else
	{
		if (totnum <= 10)
			sprintf(out, "%s%s%.1d\0", in_main, pref, pnum);
		else if (totnum <= 100)
			sprintf(out, "%s%s%.2d\0", in_main, pref, pnum);
		else
			sprintf(out, "%s%s%.3d\0", in_main, pref, pnum);
	}
}

int RedirectIOToMultipleFilesMPI(const char* fname)
{
	char filenameTMPOUT[512];

	int myrank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	AddProcessNumberToFileName(filenameTMPOUT,fname,"p_",myrank,nprocs);
	curr_stdout_fp=fopen(filenameTMPOUT,"a");
	if(curr_stdout_fp == NULL)
		return False;
#if !defined(_MSC_VER)  // TEMPORAL FIX IGOR
	*stdout=*curr_stdout_fp;
	stderr=curr_stdout_fp;
#endif
	setvbuf( stdout, NULL, _IONBF, 0 );
	ios::sync_with_stdio();
	return True;
}


void RedirectIOToConsole()
{
	return;
#if defined(_MSC_VER)
	int hConHandle;
	long lStdHandle;

	CONSOLE_SCREEN_BUFFER_INFO coninfo;
	FILE *fp;
	
	AllocConsole();
	SetConsoleTitle(TEXT("HARLEM CONSOLE"));

	int MAX_CONSOLE_LINES = 2000;

	// set the screen buffer to be big enough to let us scroll text
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE),&coninfo);
	coninfo.dwSize.Y = MAX_CONSOLE_LINES;
	SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE),coninfo.dwSize);

// redirect unbuffered STDOUT to the console

	lStdHandle = (long)GetStdHandle(STD_OUTPUT_HANDLE);
	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);

	fp = _fdopen( hConHandle, "w" );
	freopen_s(&fp, "CONOUT$", "w", stdout);

	setvbuf( stdout, NULL, _IONBF, 0 );


	lStdHandle = (long)GetStdHandle(STD_INPUT_HANDLE);
	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);

//	fp = _fdopen( hConHandle, "r" );
//	*stdin = *fp;
//	setvbuf( stdin, NULL, _IONBF, 0 );

// redirect unbuffered STDERR to the console

	lStdHandle = (long)GetStdHandle(STD_ERROR_HANDLE);
	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
	fp = _fdopen( hConHandle, "w" );
	freopen_s(&fp, "CONOUT$", "w", stderr);

	setvbuf( stderr, NULL, _IONBF, 0 );

	ios::sync_with_stdio();
#endif
}


int RedirectIOToFile(const char* fname)
{
	curr_stdout_fp=fopen(fname,"a");
	if(curr_stdout_fp == NULL)
		return False;
	*stdout=*curr_stdout_fp;
    setvbuf( stdout, NULL, _IONBF, 0 );
	ios::sync_with_stdio();
	return True;
}

int RestoreIOToConsole()
{
//	if(con_out_fp == NULL)
//		return False;
//	if(curr_stdout_fp !=NULL) fclose(curr_stdout_fp);
//		*stdout=*con_out_fp;
//		setvbuf( stdout, NULL, _IONBF, 0 );
//		ios::sync_with_stdio();
	
	return True;
}


// extern "C" DllExport 
int ErrorInMod(const char* module, const char* msg)
{
	cerr << " Error in: " << module << endl;
	cerr << msg << endl;
	return 1;
}


int ha_copy_file(const char* src, const char* tgt, const int mode )
{
	int result;
#ifdef _MSC_VER
	BOOL fail_if_exist;
	if(mode == 0) 
		fail_if_exist = FALSE;
	else
		fail_if_exist = TRUE;
	result = CopyFileA(src, tgt, fail_if_exist);
#else
	std::string cmd;
	cmd = " cp ";
	cmd += src;
	cmd += " ";
	cmd += tgt;
	system(cmd.c_str());
	result = 1;
#endif
	return result;
}

int ha_delete_file(const char* fname )
{
	int result;
#ifdef _MSC_VER
	_unlink(fname);
#else
	unlink(fname);
#endif
	result = 1;
	return result;
}


bool find_line_in_file(FILE* fp, const char* str_comp, char* cur_line, const int len, const bool rew )
{
	if(rew) rewind(fp);
	for(;;)
	{
		char* str= fgets(cur_line, len, fp);
		if(str == NULL)
			return false;
		if(!strncmp(cur_line,str_comp,strlen(str_comp)) )
			return true;
	}	
}

#if defined _MSC_VER

/*

    Implementation of POSIX directory browsing functions and types for Win32.

    Kevlin Henney (mailto:kevlin@acm.org), March 1997.

    Copyright Kevlin Henney, 1997. All rights reserved.

    Permission to use, copy, modify, and distribute this software and its
    documentation for any purpose is hereby granted without fee, provided
    that this copyright and permissions notice appear in all copies and
    derivatives, and that no charge may be made for the software and its
    documentation except to cover cost of distribution.

*/

struct DIR
{
    long                handle; /* -1 for failed rewind */
    struct _finddata_t  info;
    struct dirent       result; /* d_name null iff first time */
    char                *name;  /* NTBS */
};

DIR *opendir(const char *name)
//! Description
//! The opendir function opens the directory specified by name, which may use either / or \ as a directory separator but should not contain any wildcards. On success it associates a DIR stream with the open directory. This stream is for use in subsequent browsing operations on the directory. 
//! Returns
//! A pointer to the DIR structure for the opened directory on success, otherwise null on failure. 
//! Errors
//! ENOENT   No such directory.
//! EINVAL   Invalid argument or directory name.
//! ENOMEM   Not enough memory to perform the operation. 
{
    DIR *dir = 0;

    if(name && name[0])
    {
        size_t base_length = strlen(name);
        const char *all = /* the root directory is a special case... */
            strchr("/\\", name[base_length - 1]) ? "*" : "/*";

        if((dir = (DIR *) malloc(sizeof *dir)) != 0 &&
           (dir->name = (char *) malloc(base_length + strlen(all) + 1)) != 0)
        {
            strcat(strcpy(dir->name, name), all);

            if((dir->handle = _findfirst(dir->name, &dir->info)) != -1)
            {
                dir->result.d_name = 0;
            }
            else /* rollback */
            {
                free(dir->name);
                free(dir);
                dir = 0;
            }
        }
        else /* rollback */
        {
            free(dir);
            dir   = 0;
            errno = ENOMEM;
        }
    }
    else
    {
        errno = EINVAL;
    }

    return dir;
}

int closedir(DIR *dir)
//! Description
//! The closedir function closes the directory stream associated with dir, freeing resources as necessary and invalidating the dirpointer. 
//! Returns
//! Returns 0 on successful completion, otherwise -1. 
//! Errors
//! EBADF    Invalid directory stream. 
{
    int result = -1;

    if(dir)
    {
        if(dir->handle != -1)
        {
            result = _findclose(dir->handle);
        }

        free(dir->name);
        free(dir);
    }

    if(result == -1) /* map all errors to EBADF */
    {
        errno = EBADF;
    }

    return result;
}

struct dirent *readdir(DIR *dir)
//! Description
//! The readdir function is used to iterate through the directory stream dir. It advances it one entry at a time, details of which it returns as its result. Except for drive root directories, the caller is guaranteed that the . and .. entries will be included in the directory stream. 
//! Returns
//! Returns a pointer to the directory details on success, in which d_name is the file name of the current entry, otherwise null on error or end of stream. 
//! Errors
//! ENOENT   No more entries.
//! EBADF    Invalid directory stream.
{
    struct dirent *result = 0;

    if(dir && dir->handle != -1)
    {
        if(!dir->result.d_name || _findnext(dir->handle, &dir->info) != -1)
        {
            result         = &dir->result;
            result->d_name = dir->info.name;
        }
    }
    else
    {
        errno = EBADF;
    }

    return result;
}

void rewinddir(DIR *dir)
//! Description
//! The rewindir function can be used to reset the directory stream dir to the start. 
//! Returns
//! No error status is returned. 
//! Errors
//! EBADF    Invalid directory stream. 
{
    if(dir && dir->handle != -1)
    {
        _findclose(dir->handle);
        dir->handle = _findfirst(dir->name, &dir->info);
        dir->result.d_name = 0;
    }
    else
    {
        errno = EBADF;
    }
}

#endif
