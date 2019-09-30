from hashcode import *
from libpycommon import *
from libpyobject import *
from libpymortgl import *

import libpycommon
import libpyobject

UNKNOWN = 0

def load_mol(file, dbin=None, format=UNKNOWN):
    m = molecule_t()
    if( dbin is None ) :
        db = database_t()
    else :
        db = dbin

    libpyobject._load_mol(file, m, format)
    return m

def load_ffp(file):
    ffp = molecule_t()
    libpyobject._load_amber_ffp(file, ffp);
    return ffp;

def save_mol(file, m, format=UNKNOWN):
    libpyobject._save_mol(file, m, format)

def load_mdb(file, format=UNKNOWN):
    db = database_t()
    libpyobject._load_mdb(file, db, format)
    return db

def save_mdb(file, db, format=UNKNOWN):
    libpyobject._save_mdb(file, db, format)

def charge(obj):
    m = molecule_t()
    if type(obj) == type(m):
        return libpyobject._molcharge(obj)

    p = m.create_atom()
    if type(obj) == type(p):
        if obj.cmpid() == RESD:
	    return libpyobject._rescharge(obj)
	elif obj.cmpid() == ATOM:
	    return obj.pchg
	else:
	    print "Error: can not get charge of ", unhash( obj.cmpid() ) 
            return 0.0

    return 0.0

aliases = {"all"       : "smarts.[*]", 
           "hydrogens" : "smarts.[#1]",
	   "backbone"  : "ambmask.@C,O,N,CA"}

def mask_atom(mol, expr):
    from string import split
    (type, mask) = split(expr, '.')

    if type == "alias":
        temp = aliases[mask]
	(type, mask) = split(temp, '.')

    if type == "ambmask":
        return libpyobject._ambmsk_mask_atom(mol, mask)
    elif type == "smarts":
        return libpyobject._smarts_mask_atom(mol, mask)
    else:
        raise IndexError, "unknown type " + type

    

