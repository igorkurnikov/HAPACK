#!/usr/bin/env python

ROOT = None
PVER = None
PSUB = None
INFO = None
MTAB = "        "


def appendlog( src, dst ):
	fr = open( src, "r" )
	fw = open( dst, "a" )
	line = fr.readline()
	while len(line) > 0:
		fw.write( MTAB + MTAB + line )
		line = fr.readline()
	fr.close()
	fw.close()

def printscr( line ):
	print line

def printlog( line ):
	import os
	log = ROOT + "/src/install.log"
	os.system( "echo \"" + line + "\" >> " + log )
	
def printall( line ):
	printlog( line )
	printscr( line )

def myexec( cmd ):
	import os
	import sys
	log = ROOT + "/src/current.log"
	cmd = cmd + " >& " + log
	printall(  MTAB + cmd )
	r = os.system( cmd )
	if r != 0:
		printall( MTAB + "Failed!" )
		printall( MTAB  )
		printall( MTAB + "if it failed at wget or curl due to no internet connection, try to " )
		printall( MTAB + "download the file from other machine," )
		printall( MTAB + "copy it to the current directory and restart install.py." )
		printall( MTAB + "otherwise, check the file "+ log + " and try to resolve it" )
		sys.exit(-1)
	appendlog( log,  ROOT + "/src/install.log" )

def chdir( dir ):
	import os
	printall(  MTAB + "cd " + dir )
	os.chdir(dir)
 
def macos():
	import commands
	r = commands.getoutput( "uname" )
	return r=="Darwin"

def geturl( url, file ):
	import os
	import commands
	if macos():
		myexec( "curl " + url + " -o " + file )
	else:
		myexec( "wget " + url )

def get_pyver():
	import commands
	from string import split
	r = commands.getoutput( "python -V" )
	if r.find( "command not found" ) != -1:
		return None, None

	r = split( r )[1]
	return r[0:3], r[3:]

def get_root():
	from os import environ, getcwd
	cwd = getcwd()
	idx = cwd.find( "amber" )
	if idx != -1:
		return cwd[0:idx+7]

	if environ.has_key("AMBERHOME"):
		return environ["AMBERHOME"]

	return None

def install_clib( name, opts=None, confdir="." ):
	import os

	dict = INFO[name]

	if os.path.exists( ROOT+dict["key"] ):
		return

	printall( "Installing " + name )
	geturl( dict["url"], dict["src"] )
	myexec( "tar -zxf " + dict["src"] )
	chdir( dict["dir"] )

	config = confdir + "/configure --prefix=" + ROOT
	if not(opts is None):
		config += " " + opts;

	myexec( config )
	myexec( "make" )
	myexec( "make install" )
	chdir( ".." )

def install_leap( ):
	import os
	printall( "Installing leap" )
	build = ROOT + "/src/build"
	if os.path.exists( build ):
		os.system( "rm -rf " + build )
	os.mkdir( build )
	chdir( build )
	os.environ["AMBERHOME"]=ROOT
	myexec( ROOT + "/bin/cmake ../gleap" )
	myexec( "make" )
	myexec( "make install" )

def init():
	from os  import environ
	from sys import exit
	global ROOT, PVER, PSUB

        WELCOME = "Going to install pymort/pyleap in your system"
	print WELCOME
	print MTAB
	print "Preparation"
	print MTAB + "checking install prefix:",
	ROOT = get_root()
	if ROOT is None:
		print MTAB 
		print MTAB + "Error: I cannot locate myself. Are you sure you are running from amber source directory?"
        	print MTAB + "       Anyway, you can help me out by set environment variable AMBERHOME explicitly."
		exit(-1)
	else:
		import os
		os.system( "rm -f " + ROOT + "/src/install.log" )
		printscr( ROOT )
		printlog( WELCOME )
		printlog( MTAB + "checking install prefix: " + ROOT )

	print MTAB + "checking python version:",
	PVER,PSUB = get_pyver()
	if PVER is None:
		printall( MTAB )
		printall( MTAB + "Error: I cannot find python executable in your path" )
		exit(-2)
	else:
		printscr( "%s%s" %(PVER,PSUB) )
		printlog( MTAB + "checking python version: " + ("%s%s" %(PVER,PSUB)))


	chdir( ROOT+"/src" )
	printall( MTAB )
	printall( MTAB )

def getinf():
	global INFO
	INFO = {}
	INFO["boost"]= {"key":"/include/boost-1_34_1",
		"url":"http://superb-west.dl.sourceforge.net/sourceforge/boost/boost_1_34_1.tar.gz",
		"src":"boost_1_34_1.tar.gz",
		"dir":"boost_1_34_1"}

	INFO["cmake"]= {"key":"/bin/cmake",
		"url":"http://www.cmake.org/files/v2.6/cmake-2.6.2.tar.gz",
		"src":"cmake-2.6.2.tar.gz",
		"dir":"cmake-2.6.2"}

	INFO["python"] = {"key":"/include/python%s/Python.h"%PVER,
	  	"url":"http://www.python.org/ftp/python/%s%s/Python-%s%s.tgz" %(PVER,PSUB,PVER,PSUB),
		"src":"Python-%s%s.tgz" %(PVER,PSUB),
		"dir":"Python-%s%s" %(PVER,PSUB)}

def finish():
	printall( " " )
	printall( " " )
	printall( "Congratulations!" )
	printall( " " )
	printall( MTAB + "pymort/pyleap have been successfully installed under directory" + ROOT )
	printall( " " )
	printall( " " )
	printall( "IMPORTANT NOTES:" )
        printall( " " )
	printall( MTAB + "Please add the following directory: ")
	printall( " " )
	printall( MTAB + MTAB + ROOT + "/src/lib" )
	printall( " " )
	printall( MTAB + "to the following environmental variables: " )
	printall( " " )
	printall( MTAB + MTAB + "PYTHONPATH" )
	if macos():
		printall( MTAB + MTAB + "DYLD_LIBRARY_PATH" )
	else:
		printall( MTAB + MTAB + "LD_LIBRARY_PATH" )
	printall( " " )
	printall( "Bye!" )

init()
getinf()
install_clib( "python", "--enable-shared" )
install_clib( "boost",  "--with-libraries=python" )
install_clib( "cmake" )
install_leap( )
finish()

