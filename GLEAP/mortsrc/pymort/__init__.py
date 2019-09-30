from hashcode import *

from libpycommon import *
from libpyobject import *
from libpycapbox import ionee_i, solute_i
from libpycapbox import _addions, _addions_core
from libpycapbox import _solvatebox, _solvatebox_core
from libpycapbox import _solvatecap, _solvatecap_core
from libpycapbox import _solvateoct, _solvateoct_core
from libpycapbox import _solvateshl, _solvateshl_core
from libpyobjfun import _mortenv, get_vdwr, find_file, mdlize_mdb, mdlize_seq
from libpyformat import _load_mdb, save_top

UNKNOWN = 0

mortenv = _mortenv()

def load_mol(file, format=UNKNOWN):
	from libpyformat import _load_mol
	m = molecule_t()
	_load_mol(file, m, format)
	return m

def load_frc(file, frc=None, ftype=AMBER):
	'''
	load force field parameters.
	parameters:
		file:  file name
		frc:   the result (in type molecule_t), not None for loading frcmod
		ftype: force field type. could either be AMBER or AMOEBA
	'''
	from libpyformat import _load_frc
	if frc is None:
		frc = molecule_t()

	if ftype==AMBER:
		_load_frc(file, frc)
	elif ftype==AMOEBA:
		#_load_amoeba_frc(file, frc)
		pass
	else:
		raise RuntimeError, "Error: unknow force field type: " + unhash(ftype)

	return frc

def save_mol(m, file, format=UNKNOWN):
	from libpyformat import _save_mol
	_save_mol(file, m, format)

def load_mdb(file, format=UNKNOWN):
	db = database_t()
	_load_mdb(file, db, format)
	return db

def save_mdb(file, db, format=UNKNOWN):
	from libpyformat import _save_mdb
	_save_mdb(file, db, format)

def charge(obj):
	from libpyobject import molecule_t, morf_t
	from libpyobjfun import _molcharge, _rescharge
	from libpycapbox import ionee_i
   
	if isinstance(obj, molecule_t):
		return _molcharge(obj)

	if isinstance(obj, morf_t):
		if obj.cmpid() == RESD:
			return _rescharge(obj)
		elif obj.cmpid() == ATOM:
			return obj.pchg
		else:
			raise RuntimeError("Error: can not get charge of ", unhash(obj.cmpid()) )

	if isinstance(obj, ionee_i):
		natom = obj.natom
		pchrg = obj.getchrg()
        
		sum = 0.0
		for c in pchrg:
			sum += c

		return sum

	raise RuntimeError( "Error: cannot get charge for type " + str(type(obj)) )

aliases = {"all"       : "smarts.[*]", 
           "hydrogens" : "smarts.[#1]",
	   "backbone"  : "ambmask.@C,O,N,CA"}

def mask_atom(mol, expr):
	from string import split
	from libpyobject import _ambmsk_mask_atom
	from libpyobject import _smarts_mask_atom

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



def load_lib( name, file ):
	if not mortenv.has_key(name):
		import os
		if not os.environ.has_key("AMBERHOME"):
			raise RuntimeError( "Error: environment AMBERHOME must be set" )

		amberhome = os.environ["AMBERHOME"]
		lib = load_mdb( amberhome + "/dat/leap/lib/" + file )
		mortenv[name] = lib



def load_ion(name):
	load_lib( "_ions", "ions94.lib" );
	return mortenv["_ions"][name]

    
def load_svt(name):
	load_lib( "_solvents", "solvents.lib" )
	return mortenv["_solvents"][name]


