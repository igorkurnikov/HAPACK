#include <add.hpp>
#include <addmap.hpp>
#include <addions.hpp>
#include <bond.hpp>
#include <center.hpp>
#include <charge.hpp>
#include <copy.hpp>
#include <create.hpp>
#include <delete.hpp>
#include <desc.hpp>
#include <energy.hpp>
#include <help.hpp>
#include <impose.hpp>
#include <list.hpp>
#include <loadcrd.hpp>
#include <loadfrc.hpp>
#include <loadmol2.hpp>
#include <loadoff.hpp>
#include <loadpdb.hpp>
#include <loadprep.hpp>
#include <loadsdf.hpp>
#include <measure.hpp>
#include <merge.hpp>
#include <savemol2.hpp>
#include <saveoff.hpp>
#include <savepdb.hpp>
#include <saveprm.hpp>
#include <savesdf.hpp>
#include <set.hpp>
#include <shell.hpp>
#include <solvate.hpp>
#include <source.hpp>
#include <zmatrix.hpp>
#include <moloper.hpp>
#include <help.hpp>
#include <transform.hpp>

amber::help_command g_help_command;
amber::add_command g_add_command2;
amber::addions_command g_addions_command2;
amber::addmap_command g_addpdbatommap_command2( "addpdbatommap" );
amber::addmap_command g_addpdbresmap_command2( "addpdbresmap" );
amber::null_command g_addatomtypes_command2( "addatomtypes" );
amber::null_command g_logfile_command2( "logfile" );
amber::bond_command g_bond_command2;
amber::bondbydis_command g_bondbydis_command2;
amber::center_command g_center_command2;
amber::charge_command g_charge_command2;
amber::copy_command g_copy_command2;
amber::create_command g_createAtom_command2( "atom" );
amber::create_command g_createResd_command2(  "residue" );
amber::create_command g_createMolecule_command2( "unit" );
amber::desc_command g_desc_command2;
amber::impose_command g_impose_command2;
amber::list_command g_list_command2;    
amber::loadfrc_command g_loadamberparams_command2( "loadamberparams" );
amber::loadfrc_command g_loadamoebaparms_command2( "loadamoebaparams" );
amber::loadoff_command g_loadoff_command2;
amber::measure_command g_measure_command2;
amber::moloper_command g_fixbond_command( "fixbond" );
amber::moloper_command g_addhydr_command( "addhydr" );
amber::moloper_command g_setpchg_command( "setpchg" );
amber::moloper_command g_parmchk_command( "parmchk" );
amber::merge_command g_combine_command2( "combine" );
amber::merge_command g_sequence_command2( "sequence" );
amber::source_command g_source_command2;
amber::saveprm_command g_saveamberparm_command2( "saveamberparm" );
amber::saveprm_command g_saveamoebaprm_command2( "saveamoebaparm" );
amber::saveprm_command g_saveamberparmpol_command2( "saveamberparmpol" );
amber::saveprm_command g_savetinkerprm_command2( "savetinkerparm" );
amber::set_command g_set_command2;
amber::shell_command g_shell_command2;
amber::solvate_command g_solvatebox_command2( "box" );
amber::solvate_command g_solvateoct_command2( "oct" );
amber::solvate_command g_solvatecap_command2( "cap" );
amber::solvate_command g_solvateshl_command2( "shell" );
amber::zmatrix_command g_zmatrix_command2;
amber::loadsdf_command g_loadsdf_command2;
amber::loadpdb_command g_loadpdb_command2( "loadpdb" );
amber::loadpdb_command g_loadpdbusingseq_command2( "loadpdbusingseq" );
amber::loadmol2_command g_loadmol2_command2;
amber::savesdf_command g_savesdf_command2;
amber::savepdb_command g_savepdb_command2;
amber::saveoff_command g_saveoff_command2;
amber::savemol2_command g_savemol2_command2;
amber::energy_command g_energy_command2;
amber::loadprep_command g_loadprep_command2;
amber::transform_command translate_command2( "translate" );
amber::transform_command rotate_command2( "rotate" );
amber::transform_command transform_command2( "transform" );
