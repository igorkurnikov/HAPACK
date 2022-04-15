# __init__.py file for molset package
import os
import sys 
if( os.getenv("HARLEM_HOME") != None ): 
    os.environ["MOLSET_HOME"] = os.getenv("HARLEM_HOME")

molset_dir = os.path.dirname(__file__)
if( os.getenv("MOLSET_HOME") == None ):
    test_dirs = []
    test_dirs.append( os.path.join(molset_dir,"..") )
    test_dirs.append( os.path.join(molset_dir,"..","..","..","..","opt","interx"))
    test_dirs.append( os.path.join(molset_dir,"..","..","..","..","opt","harlem"))
    test_dirs.append( os.path.join(molset_dir,"..","..","..","..","opt","molset"))
    test_dirs.append( os.path.join(molset_dir,"..","..","..","..","share","interx"))
    test_dirs.append( os.path.join(molset_dir,"..","..","..","..","share","harlem"))
    test_dirs.append( os.path.join(molset_dir,"..","..","..","..","share","molset"))
    for d in test_dirs:
        test_res_db_dir = os.path.join(d,"residues_db")
        if( os.path.isdir( test_res_db_dir ) ):
            os.environ["MOLSET_HOME"] = d 
            os.environ["HARLEM_HOME"] = d
            break

pkg_sub_dir = "python" + str(sys.version_info[0]) + "." + str(sys.version_info[1])
pkg_dir = os.path.join(molset_dir,"..","lib",pkg_sub_dir)
sys.path.append( pkg_dir )

from molset.molsetc import *      
from molset.harlempy import start_harlem
from molset.molset_ext import *

HaAtom_FillStdAtomTypes()
HaResidue_InitStdResNames()
HaResidue_InitResSynonym()
