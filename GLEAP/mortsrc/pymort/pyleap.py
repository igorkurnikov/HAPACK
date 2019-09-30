# 
# initialization routines
#

from pymort import molecule_t, find_file, mortenv, charge

def addnmap( type, items ):
	'''
		add name map: replacement of addPdbResMap and addPdbAtomMap
	'''
	from pymort import strvec, nmap_t
	lst = strvec(len(items))
	for i in xrange(len(items)):
		lst[i] = items[i]

	if not mortenv.has_key("_namemap" ):
		mortenv["_namemap"] = nmap_t()

	nmap = mortenv["_namemap"]

	if type=="atom":
		nmap.add_atom_map( lst )
	elif type=="resd":
		nmap.add_resd_map( lst )
	else:
		raise RuntimeError, "Error: unknown type of namemap " + type

def alias( dst, src ):
	'''
		make aliases 
		example: HIS=HIE etc
	'''
	if mortenv.has_key(src):
		mortenv[dst] = mortenv[src]


def init_path():
	from os     import environ
	from pymort import mortenv
	if not environ.has_key("AMBERHOME"):
		raise RuntimeError, "environment AMBERHOME not set"
	
	root = environ["AMBERHOME"]
	path =  root + "/dat/leap/gleap:"
	path += root + "/dat/leap/parm:"
	path += root + "/dat/leap/prep:"
	path += root + "/dat/leap/lib:"
	path += root + "/dat/leap/cmd"
	mortenv.path = path

def source( file ):
	init_path()
	from string import split
	f = open( find_file(file), "r" )
	l = f.readline()
	while len(l) > 0:
		items = split( l )
		if( len(items) > 0 and items[0][0]!='#' ):
			if items[0]=="frcparm":
				loadamberparams( items[1] )
			elif items[0]=="library":
				loadoff( items[1] )
			elif items[0]=="resdmap":
				addnmap( "resd", items[1:] )
			elif items[0]=="atommap":
				addnmap( "atom", items[1:] )
			elif items[0]=="alias":
				alias( items[1], items[2] )
			else:
				raise RuntimeError, "Error: unknown keywork in initrc file " + items[1]

		l = f.readline()
#
# functions to read molecule 
#
def loadoff( file ):
	'''
		load AMBER off library into mortenv
	'''
	from pymort import _load_mdb, OFF
	file = find_file( file )
	_load_mdb( file, mortenv, OFF )

def loadpdb( file ):
	'''
		interpret pdb format, if there is residue library 
		use it to build up the molecule
	'''
	from pymort import mortenv, load_mol, mdlize_mdb, PDB
	m = load_mol( file, PDB )
	if mortenv.has_key( "ALA" ):
		mdlize_mdb( m, mortenv )
	return m

def loadsdf( file ):
	'''
		interpret mdl sdf format file
		Note: this one read only the first molecule in the file,
		to read all, use function load_mdb in pymort
	'''
	from pymort import load_mol, SDF
	return load_mol( file, SDF )

def loadmol2( file ):
	'''
		load tripos mol2 format file
		Note: this one read only the first molecule in the file,
		to read all, use function load_mdb in pymort
	'''
	from pymort import load_mol, MOL2
	return load_mol( file, MOL2 )

def loadamberparams( file ):
	''' 
		load amber force field parameters
		the force field paramter will be stored in a molecule
		_amberffp in mortenv
	'''
	from pymort import mortenv, find_file, load_frc

	if not mortenv.has_key("_amberffp"):
		mortenv["_amberffp"] = molecule_t()
		mortenv["_amberffp"].name = "amberffp"
	frc = mortenv["_amberffp"]
	file = find_file( file )
	load_frc( file, frc )

#
# functions to save molecule in various of format
#
def savesdf( m, file ):
	'''
		save molecule in MDL SDF format
	'''
	from pymort import save_mol, SDF
	save_mol( m, file, SDF )

def savemol2(m, file ):
	'''
		save molecule in tripos MOL2 format
	'''
	from pymort import save_mol, MOL2
	save_mol( m, file, MOL2 )

def saveoff( m, file ):
	'''
		save in OFF library
	'''
	from pymort import save_mol, OFF
	save_mol( m, file, OFF )

def savepdb( m, file ):
	'''
		save molecule in pdb format
	'''
	from pymort import save_mol, PDB
	save_mol( m, file, PDB )

def saveamberparm( m, topfile, xyzfile ):
	'''
		save amber toplogy file, and xyz file, 
		use the pre-loaded amber force field parameter
	'''
	from pymort import save_top, save_mol, AMBERXYZ, mortenv
	if not mortenv.has_key( "_amberffp" ):
		raise RuntimeError, "Error: no amber force field parameter was loaded"

	save_top( m, topfile, mortenv["_amberffp"] )
	save_mol( m, xyzfile, AMBERXYZ )


#
# utilty functions
#
def addions(ionee, ion, nion=-1, shlext=4.0, resolution=1.0):
	'''
		add counter ions
	'''
	from pymort import _addions, _addions_core

	if mortenv.has_key(ion):
		ion = mortenv[ion]
	else:
		raise RuntimeError, "Error: cannot find " + ion + " in the library"

	if abs(charge(ion)) < 0.05:
		raise RuntimeError, "Error: ion is neutral"

	if nion==-1:
		tchg = charge(ionee)
		if abs(tchg) < 0.05:
			raise RuntimeError, "Error: cannot neutralize the ionee since it is already neutral"

		schg = charge(ion)
		if tchg*schg > 0:
			raise RuntimeError( "Error: cannot neutralize the ionee since ion is same sign as the ionee" )

		nion = int((abs(tchg)+0.5)/abs(schg))

	if isinstance(ionee, molecule_t):
		_addions( ionee, ion, nion, shlext, resolution )
	elif isinstance(ionee, ionee_i):
		_addions_core( ionee, ion, nion, shlext, resolution )
	else:
		raise RuntimeError( "Error: cannot addions to type " + str(type(ionee)) )


def solvatebox( solute, svt, buf, clsnss=1.0 ):
	from pymort import mortenv, _solvatebox, _solvatebox_core

	if mortenv.has_key( svt ):
		svt = mortenv[svt]
	else:
		raise RuntimeError, "Error: cannot find solvent " + svt + " in the library."
	
	if isinstance(solute, molecule_t):
		return _solvatebox( solute, svt, buf, clsnss )
	elif isinstance(solute, solute_i):
		return _solvatebox_core( solute, svt, buf, clsnss )
	else:
		raise RuntimeError( "Error: cannot solvate this type: " + str(type(solute)) )

def solvateoct( solute, svt, buf, clsnss=1.0 ):
	from pymort import mortenv, _solvateoct, _solvateoct_core

	if mortenv.has_key( svt ):
		svt = mortenv[svt]
	else:
		raise RuntimeError, "Error: cannot find solvent " + svt + " in the library."
	
	if isinstance(solute, molecule_t):
		return _solvateoct( solute, svt, buf, clsnss )
	elif isinstance(solute, solute_i):
		return _solvateoct_core( solute, svt, buf, clsnss )
	else:
		raise RuntimeError( "Error: cannot solvate this type: " + str(type(solute)) )

def solvateshell( solute, svt, buf, clsnss=1.0 ):
	from pymort import mortenv, _solvateshl,  _solvateshl_core

	if mortenv.has_key( svt ):
		svt = mortenv[svt]
	else:
		raise RuntimeError, "Error: cannot find solvent " + svt + " in the library."
	
	if isinstance(solute, molecule_t):
		return _solvateshl( solute, svt, buf, clsnss )
	elif isinstance(solute, solute_i):
		return _solvateshl_core( solute, svt, buf, clsnss )
	else:
		raise RuntimeError( "Error: cannot solvate this type: " + str(type(solute)) )


