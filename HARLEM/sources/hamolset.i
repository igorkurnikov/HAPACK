%module molsetc

%begin %{
#include <mpi.h>  
#include <cstddef>
#include <Python.h>

#ifndef HA_NOGUI
#include "wx/wxprec.h" 
#include <wx/app.h>
#include <wx/docmdi.h>
#include <wx/image.h>
#include <wx/tooltip.h>  
#include <wx/filedlg.h>
#include "wx/clipbrd.h"
#include "wx/printdlg.h"

#include "wx/wx.h"
#include "wx/filename.h"  
#include <wx/splitter.h>    
#endif

%}
         
%include typemaps.i            
%include cpointer.i                      
%include stl.i                                           
%include std_iostream.i                                
%include file.i           
%include "std_string.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"

namespace std {
    %template(IntVector) vector<int>;
    %template(FloatVector) vector<float>;
	%template(DoubleVector) vector<double>;
	%template(StringVector) vector<string>;
}                    
                                                  
%{                                                           
                                                            
#include <hampi.h>                                                                   
#include "hastl.h"                 
#include "hatypes.h"                                                                    
#include "abstree.h"                     
#include "rapidxml.hpp"                                             
#include "hastring.h"                            
#include "harlemapp.h"                     
#include "hampi.h"   
#include "haconst.h"
#include "haconsts.h"
#include "haxml.h" 
#include "tinyxml.h"
#include "halinalg.h"        
#include "haatom.h"  
#include "haatgroup.h"
%}

%template(AtomGroupList) std::list<AtomGroup>;
%template(HaAtomVector) std::vector<HaAtom*>;
%template(CrdSnapshotVector) std::vector<CrdSnapshot*>;
%template(AtomIntMap) std::map<HaAtom*, int, less<HaAtom*> >;

%{
#include "hamolecule.h"          
#include "atom_mapping.h" 
#include "hamolset.h"  
#include "moleditor.h"   
#include "haresdb.h"
#include "qc_params.h"
#include "haqchem.h"
#include "hamultipole.h" 
#include "hapseudopot.h"
#include "hagaussian.h"
#include "hazindo.h"
#include "haatbasdb.h"  
#include "halocorb.h"
#include "hacoord.h"
#include "haintcrd.h"
#include "rigidbodycoord.h"
#include "trajanal.h"
#include "haenefunc.h"
#include "haintermol.h"
#include "etcoupl.h"
#include "hamatdb.h"
#include "hatypes.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_params.h"
#include "hamolmech.h" 
#include "mm_traj_anal.h"
#include "mm_driver_amber.h"
#include "mm_driver_tinker.h"
#include "mm_force_field.h"
#include "electrostmod.h"
#include "elmod.h"
#include "nuclacidmod.h"
#include "stmmod.h"
#include "canvas3d.h"
#include "hamolview.h"
#include "haempirical.h"
#include "hasurface.h"
#include "apbsmod.h"
#include "hatests.h"
#include "moleditor.h" 
#include "protonredox.h"
#ifndef HA_NOGUI
#include "dialogs_wx_2.h"
void StartHaMainFrameWX();
#endif
void StartHarlemApp();
%}

%pointer_class(int,intp);
%pointer_class(double,doublep);
%pointer_cast(int*,void*,intp_to_voidp)
%pointer_cast(void*,int*,voidp_to_intp)

%include "abstree.h"
%include "haconst.h"

%include "haio.h"
%include "hastl.h"
%include "hastring.h"  
#include "hatypes.h"
%include "haconsts.h"
%include "tinyxml.h"
%include "haxml.h"
%include "haobject.h"
%include "halinalg.h"
%include "command.h"
%include "abstree.h"
%include "vec3d.h"  
%include "haatom.h"
%include "haatgroup.h"
%include "habond.h"
%include "object3d.h"
%include "hamolecule.h"
%include "atom_mapping.h"
%include "hamolset.h"
%include "harlemapp.h"
%include "hampi.h"

%include "canvas3d.h"
%include "hamolview.h"
%include "hasurface.h"

%include "haatombasis.h"
%include "halocorb.h"
%include "hamultipole.h"
%include "hapseudopot.h" 
%include "hacompmod.h"
%include "qc_params.h"
%include "haqchem.h"    
%include "etcoupl.h"
%include "electrostmod.h"
%include "elmod.h"
%include "apbsmod.h"
%include "hagaussian.h"
%include "hazindo.h"
%include "hacoord.h" 
%include "haintcrd.h" 
%include "rigidbodycoord.h"
%include "trajanal.h"
%include "haenefunc.h"
%include "hasimulator.h"
%include "haintermol.h" 
%include "haempirical.h"
%include "hatypes.h" 
%include "mm_elements.h"
%include "mm_model.h"
%include "mm_params.h"
%include "mm_traj_anal.h"
%include "mm_force_field.h"
%include "mm_driver_amber.h"
%include "hamolmech.h"
%include "nuclacidmod.h"

%include "haresdb.h"
%include "hamatdb.h"
%include "hatests.h"
%include "moleditor.h" 
%include "protonredox.h"

#ifndef HA_NOGUI
%include "dialogs_wx_2.h"

void StartHaMainFrameWX();
#endif

void StartHarlemApp();
