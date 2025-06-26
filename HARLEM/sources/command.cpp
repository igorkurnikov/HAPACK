/*! \file command.cpp

   \author Igor Kurnikov
   \date 1997-2003

   based on 
   command.c
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
*/
#define COMMAND_CPP

#include <mpi.h>

#include "haconst.h" 
#include "Python.h"

#include <boost/algorithm/string.hpp>

#include "haio.h"

#if !defined(_MSC_VER)
#include <pwd.h>
#endif

#include "harlemapp.h"
#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "hamatdb.h"
#include "etcoupl.h"

#include "command.h"
#include "tokens.h"
#include "abstree.h"
#include "hamolview.h"

#include "haqchem.h"
#include "hagaussian.h"
#include "hadalton.h"




#define IsIdentChar(x)  ((isalnum(x))||((x)=='_')||((x)=='$'))

int CmdParser::max_history_save=1000;

StrIntMap CmdParser::Keywords;

CmdParser::CmdParser()
{
  CurToken=0;
  TokenValue=0;
  cursor_pos = 0;
  cur_history_cmd=-1;
  TokenValueFloat=0.0;
  TokenIdent = "";
  str_start_pos=0;
}

CmdParser::~CmdParser()
{
	
}

int CmdParser::InitKeywords()
{
	RegisterKeyword( "X", XTok );
	RegisterKeyword( "Y", YTok );
	RegisterKeyword( "Z" , ZTok );
	RegisterKeyword( "AT", ATTok );
	RegisterKeyword( "CG", CGTok );
	RegisterKeyword( "DB", DbTok );
	RegisterKeyword( "ET", ETTok );
	RegisterKeyword( "MO", MOTok );
	RegisterKeyword( "ON", TrueTok );
	RegisterKeyword( "OR", OrTok );
	RegisterKeyword( "ALL", AllTok );
	RegisterKeyword( "AND", AndTok );
	RegisterKeyword( "BMP", BMPTok );
	RegisterKeyword( "CEX", CEXTok );
	RegisterKeyword( "CIF", CIFTok );
	RegisterKeyword( "CPK", CPKTok );
	RegisterKeyword( "DNA", DNATok );
	RegisterKeyword( "GIF", GIFTok );
	RegisterKeyword( "ION", IonTok );
	RegisterKeyword( "KEY", KeyTok );
	RegisterKeyword( "MDL", MDLTok );
	RegisterKeyword( "NOT", NotTok );
	RegisterKeyword( "OFF" , FalseTok );
	RegisterKeyword( "PDB" , PDBTok );
	RegisterKeyword( "PPM" , PPMTok );
	RegisterKeyword( "RED" , RedTok );
	RegisterKeyword( "RGB" , IRISTok);
	RegisterKeyword( "RNA" , RNATok );
	RegisterKeyword( "RUN" , RunTok );
	RegisterKeyword( "RWF" , RWFTok );
	RegisterKeyword( "SET" , SetTok );
	RegisterKeyword( "SUN" , SUNTok );
	RegisterKeyword( "VDW" , VDWTok );
	RegisterKeyword( "XYZ" , XYZTok );
	RegisterKeyword( "ZAP" , ZapTok );
	RegisterKeyword( "ATOM" , AtomTok );
	RegisterKeyword( "AXES" , AxesTok );
	RegisterKeyword( "BLUE" , BlueTok );
	RegisterKeyword( "BOND" , BondTok );
	RegisterKeyword( "CYAN" , CyanTok );
	RegisterKeyword( "DASH" , DashTok );
	RegisterKeyword( "DOTS" , DotsTok );
	RegisterKeyword( "DUMP" , DumpTok );
	RegisterKeyword( "ECHO" , EchoTok );
	RegisterKeyword("EDIT", EditTok);	
	RegisterKeyword("EPSF", EPSFTok );
	RegisterKeyword("EXIT", ExitTok );
	RegisterKeyword("FDAT", FDATTok );
	RegisterKeyword("HALF", HalfTok );
	RegisterKeyword("HELP", HelpTok );
	RegisterKeyword("IDX2", Idx2Tok );
	RegisterKeyword("IDX4", Idx4Tok  );
	RegisterKeyword("INFO", InfoTok );
	RegisterKeyword("INIT", InitTok );
	RegisterKeyword("IONS", IonTok  );
	RegisterKeyword("IRIS", IRISTok );
    RegisterKeyword("JPEG", JPEGTok );
	RegisterKeyword("LIST", ListTok );
	RegisterKeyword("LOAD", LoadTok );
	RegisterKeyword("MMDB", MMDBTok  );
	RegisterKeyword("MODE", ModeTok );
	RegisterKeyword("MOL2", Mol2Tok );
	RegisterKeyword("MONO", MonoTok );
	RegisterKeyword("MULT", MultTok );
	RegisterKeyword("NAME", NameTok );
	RegisterKeyword("NONE", NoneTok );
	RegisterKeyword("OPEN", OpenTok );
	RegisterKeyword("OPER", OperTok );
	RegisterKeyword("PICT", PICTTok );
	RegisterKeyword("QUIT", QuitTok );
	RegisterKeyword("READ", ReadTok );
//	RegisterKeyword("SAVE", SaveTok);
	RegisterKeyword("SHOW", ShowTok );
	RegisterKeyword("SLAB", SlabTok );
//	RegisterKeyword("TEST", TestTok );
	RegisterKeyword("TRUE", TrueTok );
	RegisterKeyword("TURN", TurnTok );
	RegisterKeyword("TYPE", TypeTok );
	RegisterKeyword("USER", UserTok );
	RegisterKeyword("VRML", VRMLTok );
	RegisterKeyword("WAIT", WaitTok );
	RegisterKeyword("ZOOM", ZoomTok);
	RegisterKeyword("ALPHA", AlphaTok );
	RegisterKeyword("AMINO", AminoTok );
	RegisterKeyword("ANGLE", AngleTok );
	RegisterKeyword("ATOMS", AtomTok  );
	RegisterKeyword("BASIC", BasicTok );
	RegisterKeyword("BASIS", BasisTok );
	RegisterKeyword("BLACK", BlackTok );
	RegisterKeyword("BONDS", BondTok  );
	RegisterKeyword("CHAIN", ChainTok );
	RegisterKeyword("COLOR", ColourTok );
	RegisterKeyword("COORD", CoordTok );
	RegisterKeyword("DONOR", DonorTok);
	RegisterKeyword("FALSE", FalseTok);
	RegisterKeyword("GREEN", GreenTok );
	RegisterKeyword("GROUP", GroupTok );
	RegisterKeyword("HBOND", HBondTok );
	RegisterKeyword("HELIX", HelixTok  );
	RegisterKeyword("IDENT", IdentifyTok );
	RegisterKeyword("LABEL", LabelTok );
	RegisterKeyword("LARGE", LargeTok  );
	RegisterKeyword("MENUS", MenusTok  );
	RegisterKeyword("MODEL", ModelTok );
	RegisterKeyword("MOPAC", MOPACTok );
	RegisterKeyword("MOUSE", MouseTok  );
	RegisterKeyword("PAUSE", WaitTok  );
	RegisterKeyword("POLAR", PolarTok );
	RegisterKeyword("PRINT", PrintTok );
	RegisterKeyword("QCHEM", QChemTok );
	RegisterKeyword("RENUM", RenumTok );
	RegisterKeyword("RESET", ResetTok  );
	RegisterKeyword("RESNO", ResNoTok );
	RegisterKeyword("SHEET", SheetTok );
	RegisterKeyword("SHELX", SHELXTok );
	RegisterKeyword("SMALL", SmallTok  );
	RegisterKeyword("SOLID", SolidTok );
	RegisterKeyword("SYBYL", SybylTok );
	RegisterKeyword("TRACE", TraceTok );
	RegisterKeyword("TURNS", TurnTok );
	RegisterKeyword("WATER", WaterTok );
	RegisterKeyword("WHITE", WhiteTok );
	RegisterKeyword("ACIDIC", AcidicTok );
	RegisterKeyword("ANGLES", AngleTok );
	RegisterKeyword("ATOMNO", AtomNoTok );
	RegisterKeyword("BIOSYM", BiosymTok );
	RegisterKeyword("BONDED", BondedTok );
	RegisterKeyword("BURIED", BuriedTok );
	RegisterKeyword("CENTER", CentreTok );
	RegisterKeyword("CENTRE", CentreTok );
	RegisterKeyword("CHARGE", ChargeTok );
	RegisterKeyword("CHARMM", CharmmTok );
	RegisterKeyword("CHOOSE", ChooseTok );
	RegisterKeyword("COLORS", ColourTok );
	RegisterKeyword("COLOUR", ColourTok );
	RegisterKeyword("CREATE", CreateTok );
	RegisterKeyword("CYCLIC", CyclicTok );
	RegisterKeyword("DALTON", DaltonTok );
	RegisterKeyword("DASHES", DashTok );
	RegisterKeyword("DEFINE", DefineTok );
	RegisterKeyword("DELETE", DeleteTok );
	RegisterKeyword("ELEMNO", ElemNoTok );
	RegisterKeyword("ENERGY", EnergyTok );
	RegisterKeyword("GROUPS", GroupsTok );
	RegisterKeyword("HARLEM", HarlemTok );
	RegisterKeyword("HBONDS", HBondTok );
	RegisterKeyword("HETERO", HeteroTok );
	RegisterKeyword("HOLLOW", HollowTok );
	RegisterKeyword("INLINE", InLineTok );
	RegisterKeyword("LABELS", LabelTok );
	RegisterKeyword("LIGAND", LigandTok );
	RegisterKeyword("LOCORB", LocOrbTok );
	RegisterKeyword("MEDIUM", MediumTok );
	RegisterKeyword("MODULE", ModuleTok );
	RegisterKeyword("MONOPS", MonoPSTok );
	RegisterKeyword("NMRPDB", NMRPDBTok );
	RegisterKeyword("NORMAL", NormalTok );
	RegisterKeyword("ORANGE", OrangeTok );
	RegisterKeyword("OUTPUT", OutputTok );
	RegisterKeyword("POVRAY", POVRayTok );
	RegisterKeyword("PURINE", PurineTok );
	RegisterKeyword("PURPLE", PurpleTok );
	RegisterKeyword("QUANTA", QuantaTok );
	RegisterKeyword("RADIUS", RadiusTok );
	RegisterKeyword("RASMOL", RasMolTok );
	RegisterKeyword("RASWIN", RasMolTok );
	RegisterKeyword("REJECT", RejectTok );
	RegisterKeyword("RESIZE", ResizeTok );
	RegisterKeyword("RIBBON", RibbonTok );
	RegisterKeyword("ROTATE", RotateTok );
	RegisterKeyword("SCRIPT", ScriptTok );
	RegisterKeyword("SELECT", SelectTok );
	RegisterKeyword("SHADOW", ShadowTok );
	RegisterKeyword("SHEETS", SheetTok );
	RegisterKeyword("SOURCE", SourceTok );
	RegisterKeyword("SSBOND", SSBondTok );
	RegisterKeyword("STEREO", StereoTok );
	RegisterKeyword("SUNRLE", SUNRLETok );
	RegisterKeyword("UPDATE", UpdateTok );
	RegisterKeyword("VECTPS", VectPSTok );
	RegisterKeyword("VIOLET", VioletTok );
	RegisterKeyword("WATERS", WaterTok );
	RegisterKeyword("WITHIN", WithinTok );
	RegisterKeyword("YELLOW", YellowTok  );
	RegisterKeyword("ACYCLIC", AcyclicTok );
	RegisterKeyword("ALCHEMY", AlchemyTok );
	RegisterKeyword("AMBIENT", AmbientTok );
	RegisterKeyword("CHEMGRP", ChemGrpTok );
	RegisterKeyword("CARTOON", CartoonTok );
	RegisterKeyword("CHARGED", ChargedTok );
	RegisterKeyword("CHARGES", ChargeTok );
	RegisterKeyword("COLOURS", ColourTok );
	RegisterKeyword("CONNECT", ConnectTok );
	RegisterKeyword("CONSOLE", ConsoleTok );
	RegisterKeyword("CYSTINE", CystineTok );
	RegisterKeyword("DISPLAY", DisplayTok );
	RegisterKeyword("HELICES", HelixTok  );
	RegisterKeyword("INPFILE", InpFileTok );
	RegisterKeyword("INSIGHT", InsightTok );
	RegisterKeyword("LIGANDS", LigandTok );
	RegisterKeyword("MAGENTA", MagentaTok );
	RegisterKeyword("MOLINFO", MolinfoTok  );
	RegisterKeyword("MOLMECH", MolMechTok );
	RegisterKeyword("MONITOR", MonitorTok );
	RegisterKeyword("NEUTRAL", NeutralTok );
	RegisterKeyword("NUCLEIC", NucleicTok );
	RegisterKeyword("OVERLAP", OverlapTok );
	RegisterKeyword("PICKING", PickingTok );
	RegisterKeyword("PROTEIN", ProteinTok );
	RegisterKeyword("PURINES", PurineTok  );
	RegisterKeyword("REFRESH", RefreshTok );
	RegisterKeyword("RESIDUE", ResidueTok );
	RegisterKeyword("RIBBON1", Ribbon1Tok );
	RegisterKeyword("RIBBON2", Ribbon2Tok );
	RegisterKeyword("RIBBONS", RibbonTok );
	RegisterKeyword("SECTION", SectionTok );
	RegisterKeyword("SHADOWS", ShadowTok );
	RegisterKeyword("SHAPELY", ShapelyTok  );
	RegisterKeyword("SOLVENT", SolventTok );
	RegisterKeyword("SSBONDS", SSBondTok );
	RegisterKeyword("STRANDS", StrandsTok );
	RegisterKeyword("SURFACE", SurfaceTok );
	RegisterKeyword("TORSION", TorsionTok );
	RegisterKeyword("WAVEFUN", WaveFunTok );
	RegisterKeyword("ACCEPTOR", AcceptorTok );
	RegisterKeyword("AROMATIC", AromaticTok );
	RegisterKeyword("BACKBONE", BackboneTok );
	RegisterKeyword("BACKFADE", BackFadeTok );
	RegisterKeyword("BESTPATH", BestPathTok );
	RegisterKeyword("BONDMODE", BondModeTok );
	RegisterKeyword("BOUNDBOX", BoundBoxTok );
	RegisterKeyword("CARTOONS", CartoonTok  );
	RegisterKeyword("COUPLMAP", CouplMapTok );
	RegisterKeyword("DEPTHCUE", DepthCueTok );
	RegisterKeyword("DISTANCE", DistanceTok );
	RegisterKeyword("ELSTATIC", ElStaticTok );
	RegisterKeyword("FONTSIZE", FontSizeTok );
	RegisterKeyword("GAUSSIAN", GaussianTok );
	RegisterKeyword("GRPCONTR", GrpContrTok );
	RegisterKeyword("HARDWARE", HardwareTok );
	RegisterKeyword("HYDROGEN", HydrogenTok);
	RegisterKeyword("IDENTIFY", IdentifyTok );
	RegisterKeyword("INTERMOL", InterMolTok );
	RegisterKeyword("MOLECULE", MoleculeTok );
	RegisterKeyword("MONITORS", MonitorTok );
	RegisterKeyword("NEGATIVE", AcidicTok );
	RegisterKeyword("PATHWAYS", PathwaysTok );
	RegisterKeyword("POSITIVE", BasicTok  );
	RegisterKeyword("RENUMBER", RenumTok );
	RegisterKeyword("RESTRICT", RestrictTok  );
	RegisterKeyword("RIBBONS1", Ribbon1Tok );
	RegisterKeyword("RIBBONS2", Ribbon2Tok );
	RegisterKeyword("ROTANGLE", RotAngleTok );
	RegisterKeyword("SELECTED", SelectedTok );
	RegisterKeyword("SEQUENCE", SequenceTok );
	RegisterKeyword("SLABMODE", SlabModeTok  );
	RegisterKeyword("SOLVENTS", SolventTok  );
	RegisterKeyword("SPECULAR", SpecularTok );
	RegisterKeyword("SYMMETRY", SymmetryTok );
	RegisterKeyword("TORSIONS", TorsionTok );
	RegisterKeyword("UNITCELL", UnitCellTok );
	RegisterKeyword( "ALIPHATIC", AliphaticTok);
	RegisterKeyword( "CALCULATE", CalculateTok );
	RegisterKeyword( "CLIPBOARD", ClipboardTok);
	RegisterKeyword( "DISTANCES", DistanceTok);
	RegisterKeyword( "GREENBLUE", GreenblueTok);
	RegisterKeyword( "HOURGLASS", HourGlassTok);
	RegisterKeyword( "MAINCHAIN", MainChainTok);
	RegisterKeyword( "MOLSCRIPT", MolScriptTok);
	RegisterKeyword( "MOUSEMODE", MouseTok);
	RegisterKeyword( "POTENTIAL", PotentialTok);
	RegisterKeyword( "REDORANGE", RedorangeTok);
	RegisterKeyword( "SIDECHAIN", SidechainTok);
	RegisterKeyword( "SPACEFILL", SpacefillTok);
	RegisterKeyword( "SPECPOWER", SpecPowerTok );
	RegisterKeyword( "STRUCTURE", StructureTok);
	RegisterKeyword( "TRANSLATE", TranslateTok);
	RegisterKeyword( "WIREFRAME", WireframeTok);
	RegisterKeyword( "BACKGROUND", BackgroundTok );
	RegisterKeyword( "FOREGROUND", ForegroundTok);
	RegisterKeyword( "MACROMODEL", MacroModelTok );
	RegisterKeyword( "MONOCHROME", MonoTok );
	RegisterKeyword( "OVERLAPMOL", OverlapMolTok);
	RegisterKeyword( "PYRIMIDINE", PyrimidineTok);
	RegisterKeyword( "BOUNDINGBOX", BoundBoxTok);
	RegisterKeyword( "HYDROPHOBIC", HydrophobicTok);
	RegisterKeyword( "INFORMATION", InfoTok);
	RegisterKeyword( "PYRIMIDINES", PyrimidineTok);
	RegisterKeyword( "TEMPERATURE", TemperatureTok);
	RegisterKeyword( "TRANSPARENT", TransparentTok);
	RegisterKeyword( "SAVEAMBERPARM", SaveAmberParmTok);
	RegisterKeyword( "ALIGNOVERLAPMOL", AlignOverlapMolTok);
	RegisterKeyword( "LOADAMBERRESTART", LoadAmberRestartTok);
	RegisterKeyword( "COMBINEBONDINGATOMS", CombineBondingAtomsTok);
	
	return TRUE;
}

int CmdParser::RegisterKeyword(std::string keyw, const int itok)
{
	if(Keywords.find(keyw) != Keywords.end() )
	{
		char buf[256];
		sprintf(buf," Keyword %s already registered ",keyw.c_str());
		ErrorInMod("CmdParser::RegisterKeyword()", buf);
		return FALSE;
	}
	Keywords[keyw]=itok;
	return TRUE;
}

int CmdParser::SetCmdLine(const std::string& cmd_line)
{
	CurLine = cmd_line;
	boost::trim(CurLine);

	ResetCursorPosition();
	cmd_history.push_front(CurLine);
	if(cmd_history.size() > max_history_save)
	{
		cmd_history.pop_back();
	}
	cur_history_cmd=0;
	return TRUE;
}

const char* CmdParser::GetCmdLine()
{
	return CurLine.c_str();
}

const char* CmdParser::GetStartPosSubstr()
{
	if( str_start_pos < CurLine.size())
	{
		return &CurLine[str_start_pos];
	}
	return CurLine.c_str();
}

std::string CmdParser::RollHistory(int step)
{
	cur_history_cmd += step;
	if(cur_history_cmd < 0 || (cmd_history.size() == 0))
	{
		cur_history_cmd = -1;
		return "";
	}

	if( cur_history_cmd > ( cmd_history.size() - 1) )
	{
		cur_history_cmd= cmd_history.size();
		return "";
	}

	return cmd_history[cur_history_cmd];

}

void CmdParser::ResetCursorPosition()
{
	cursor_pos = 0;
	str_start_pos = 0;
}

void CmdParser::CommandError(const char* error )
{
    if( cursor_pos >= CurLine.size() || CurLine[cursor_pos] == 0 )
    {
		PrintLog("        \n");
		int i;
        for( i = 0; i < str_start_pos; i++ )
            PrintLog("%c",' ');
        PrintLog("^\n");
    }
		
    if( error )
    {
		PrintLog(error);
        PrintLog("!\n");
    }
    CurToken = 0;
}

int MolSet::FetchFile(int format, std::string file_name, const AtomLoadOptions& opt_par )
{
	AtomLoadOptions opt(opt_par);

    int done;
	std::string fname = file_name;
	boost::trim(fname);

	if(GetNMol() == 0)
	{
		std::string mol_name= harlem::GetPrefixFromFullName(fname);
		SetName(mol_name.c_str());
	}
	
	switch( format )
	{	
	case(FormatHarlem):      done = LoadHarlemFile(fname,    opt); break;
	case(FormatAmberPrep):   done = LoadAmberPrepFile(fname, opt); break;
	case(FormatAmberTop):    done = LoadAmberTopFile (fname, opt); break;
	case(FormatAmberOff):    done = LoadAmberOffFile(fname,  opt); break;
	case(FormatRWF):         done = LoadRWFMolecule(fname,   opt); break;
	case(FormatPDB):         done = LoadPDBFile(fname,       opt); break;
	case(FormatNMRPDB):      opt.set_i("NMRPDB", 1);  done = LoadPDBFile(fname, opt); break;
	case(FormatMol2):        done = LoadMol2File(fname,      opt); break;
	case(FormatMDL):         done = LoadMDLFile(fname,       opt); break;
	case(FormatXYZ):         done = LoadXYZFile(fname,       opt); break;
	case(FormatHIN):         done = LoadHINFile(fname,       opt); break;
	case(FormatNRG):         done = LoadNRGFile(fname,       opt); break;

	default:                 done = False;
	}
		
    if( !done )
    {
		return( False );
    }
	
	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
        (*mol_itr)->DescribeMolecule();
	}
		
	HaMolView* pView= GetActiveMolView();
	if(!pView)
		return True;

	if( GetNHBonds() > 0 )
		pView->SetHBondStatus(True,0);
	

    pView->ReDrawFlag |= (RFInitial);
    pView->InitialTransform();
	    
    return( done );
}


int CmdParser::LookUpKeyword()
{
	map<std::string,int, less<std::string> >::iterator kitr;

	boost::to_upper(TokenIdent);

	kitr = Keywords.find(TokenIdent);
	if( kitr == Keywords.end())
		return( IdentTok);

	return( (*kitr).second);
	
}

int CmdParser::FetchToken()
{
    char ch;

	int line_size = CurLine.size();

    CurToken = 0;
    while( True )
    {
	// Scroll the pointer to the first non-space character
    // or return zero if end of string or # sign is encountered

		if(cursor_pos >= line_size ) return 0;

		ch = CurLine[cursor_pos];
		cursor_pos++;

		if( !ch || (ch=='#') ) return(0);
		if( isspace(ch) ) continue;
		
		str_start_pos = cursor_pos-1;
		if( isalpha(ch) )
		{
			TokenIdent = "";
			TokenIdent += toupper(ch);
			while( IsIdentChar( CurLine[cursor_pos] ) )
			{
				ch = CurLine[cursor_pos];
				TokenIdent += toupper(ch);
				cursor_pos++;
				if( cursor_pos == line_size ) break;
			}
			return( CurToken = LookUpKeyword() );
		}
		else if( isdigit(ch) || ch == '-' || ch == '.')
		{
			double scale = 0.0;
			int isign = 1;
			CurToken = NumberTok;
			TokenValue = 0;

			if(ch == '-')
			{
				TokenValue=0;
				isign=-1;
			}

			else if(ch == '.')
			{
				cursor_pos--;
//				CurToken = FloatTok;
			}
			else
			{
				TokenValue = ch-'0';
			}

			for(;;)
			{
				if( CurLine[cursor_pos] == '.' )
				{
					CurToken = FloatTok;
					scale = 0.1;

					TokenValueFloat = TokenValue;
					cursor_pos++;
				}
				else if(isdigit(CurLine[cursor_pos]))
				{
					int val = CurLine[cursor_pos]-'0';
					if(CurToken == NumberTok)
					{
						TokenValue = 10*TokenValue + val;
					}
					else if(CurToken == FloatTok)
					{
						TokenValueFloat+= scale* (double)val;
						scale*= 0.1;
					}
					cursor_pos++;
					if( cursor_pos == line_size ) break;
				}
				else
				{
					break;
				}
			}

			if( CurToken == NumberTok )
			{
				TokenValue *= isign;
				TokenValueFloat = TokenValue;
			}
			else if( CurToken == FloatTok )
			{
				TokenValueFloat*= (double) isign;
			}
			
			return( CurToken );
		}
		else if( (ch=='\'') || (ch=='\"') || (ch=='`') )
		{
			TokenIdent = "";
			while( CurLine[cursor_pos] != 0 && (CurLine[cursor_pos] != ch) && cursor_pos < CurLine.size() )
			{
				TokenIdent += CurLine[cursor_pos];
				cursor_pos++;
			}
			if( cursor_pos >= CurLine.size() || CurLine[cursor_pos] == 0)
			{
				CommandError("String constant unterminated");
				return 0;
			}
			cursor_pos++;
			return( CurToken = StringTok );
		}
		else if( ispunct(ch) )
		{
			return( CurToken = ch );
		}
    }
}


int CmdParser::NextIf(int token, const char* error )
{
    if( FetchToken()!=token )
    {
		PrintLog(error);
		PrintLog("\n");
        return( True );
    }
	else
		return( False );
}


int CmdParser::ParseColour(int& RVal, int& GVal, int& BVal)
{
    switch( CurToken )
    {
	    case(BlueTok):        RVal=0;   GVal=0;   BVal=255; break;
        case(BlackTok):       RVal=0;   GVal=0;   BVal=0;   break;
        case(CyanTok):        RVal=0;   GVal=255; BVal=255; break;
        case(GreenTok):       RVal=0;   GVal=255; BVal=0;   break;
        case(GreenblueTok):   RVal=46;  GVal=139; BVal=87;  break;
        case(MagentaTok):     RVal=255; GVal=0;   BVal=255; break;
        case(OrangeTok):      RVal=255; GVal=165; BVal=0;   break;
        case(PurpleTok):      RVal=160; GVal=32;  BVal=240; break;
        case(RedTok):         RVal=255; GVal=0;   BVal=0;   break;
        case(RedorangeTok):   RVal=255; GVal=69;  BVal=0;   break;
        case(VioletTok):      RVal=238; GVal=130; BVal=238; break;
        case(WhiteTok):       RVal=255; GVal=255; BVal=255; break;
        case(YellowTok):      RVal=255; GVal=255; BVal=0;   break;

        case('['):    RVal = GVal = BVal = 0;

                      if( NextIf(NumberTok,"Integer value expected") ) { return(False);
                      } else if( TokenValue>255 )
                      {   PrintLog("Parameter value too large"); return( False );
                      } else RVal = (int)TokenValue;

                      if( NextIf(',',"Comma separator missing") ) return(False);
                      if( NextIf(NumberTok,"Integer value expected") ) { return(False);
                      } else if( TokenValue>255 )
                      {   PrintLog("Parameter value too large"); return( False );
                      } else GVal = (int)TokenValue;

                      if( NextIf(',',"Comma separator missing") ) return(False);
                      if( NextIf(NumberTok,"Integer value expected") ) { return(False);
                      } else if( TokenValue>255 )
                      {   PrintLog("Parameter value too large"); return( False );
                      } else BVal = (int)TokenValue;

                      return( !NextIf(']',"Close bracket ']' expected") );

        default:  return(False);
    }
    return( True );
}



AtomExpr* CmdParser::ParseRange(int neg )
{
    AtomExpr *tmp1,*tmp2;
    char ch;
	
    tmp1 = new AtomExpr();
    tmp1->type = OpLftProp|OpRgtVal;
    tmp1->rgt.val = neg? -(int)TokenValue : (int)TokenValue;
    tmp1->lft.val = PropResId;

	int line_size = CurLine.size();
	
    if( cursor_pos < line_size && CurLine[cursor_pos] == '-' )
    {
		cursor_pos++;
        neg = (CurLine[cursor_pos] == '-');
        if( neg ) cursor_pos++;
        FetchToken();
		
        if( CurToken != NumberTok )
        {
			PrintLog("Integer value expected\n");;
			delete tmp1;
            return( NULL );
        }
		
        tmp1->type |= OpMoreEq;
        tmp2 = new AtomExpr();
        tmp2->rgt.ptr = tmp1;
        tmp2->type = OpAnd;
		
        tmp1 = new AtomExpr();
        tmp1->type = OpLftProp|OpRgtVal|OpLessEq;
        tmp1->rgt.val = neg? -(int)TokenValue : (int)TokenValue;
        tmp1->lft.val = PropResId;
        tmp2->lft.ptr = tmp1;
        tmp1 = tmp2;
    }
	else
	{
		tmp1->type |= OpEqual;
	}
	
    if( cursor_pos < line_size && CurLine[cursor_pos] == ':' )
	{
        cursor_pos++;
	}
	
	if( cursor_pos < line_size )
	{
		ch = CurLine[cursor_pos];
	}
	else
	{
		ch = 0;
	}
    if( isalnum(ch) )
    {
		ch = toupper(ch);
        cursor_pos++;
		
        tmp2 = new AtomExpr();
        tmp2->type = OpAnd;
        tmp2->rgt.ptr = tmp1;
		
        tmp1 = new AtomExpr();
        tmp1->type = OpEqual | OpLftProp | OpRgtVal;
        tmp1->lft.val = PropChain;
        tmp1->rgt.val = ch;
		
        tmp2->lft.ptr = tmp1;
        tmp1 = tmp2;
    }
	else if( (ch=='?') || (ch=='%') || (ch=='*') )
	{
        cursor_pos++;
	}
	
    FetchToken();
    return( tmp1 );
}


AtomExpr* CmdParser::ParseExpression(int level, MolSet* pmset )
{
    AtomExpr *tmp1,*tmp2;
    int done, pred;
    int neg;
	double dtmp;

	vector<HaMolecule*>::iterator mol_itr;
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	
    switch( level )
    {
	case(0): /* Disjunctions */
		tmp1 = ParseExpression(1,pmset);
		while( (CurToken==OrTok) || (CurToken=='|') || (CurToken==',') )
		{
			if( CurToken=='|' )
			{
				if( FetchToken()=='|' ) FetchToken();
			}
			else
			{
				FetchToken();
			}
			
			tmp2 = new AtomExpr();
			tmp2->type = OpOr;
			tmp2->lft.ptr = tmp1;
			tmp2->rgt.ptr = NULL;
			if( !(tmp1=ParseExpression(1,pmset)) )
			{
				delete tmp2;
				return( tmp1 );
			}
			tmp2->rgt.ptr = tmp1;
			tmp1 = tmp2;
		}
		return( tmp1 );
		
	case(1): /* Conjunctions */
		tmp1 = ParseExpression(2,pmset);
		while( (CurToken==AndTok) || (CurToken=='&') )
		{
			if( CurToken=='&' )
			{
				if( FetchToken()=='&' ) FetchToken();
			}
			else
			{
				FetchToken();
			}
			
			tmp2 = new AtomExpr;
			tmp2->type = OpAnd;
			tmp2->lft.ptr = tmp1;
			tmp2->rgt.ptr = NULL;
			if( !(tmp1=ParseExpression(2,pmset)) )
			{
				delete tmp2;
				return( tmp1 );
			}
			tmp2->rgt.ptr = tmp1;
			tmp1 = tmp2;
		}
		return( tmp1 );
		
	case(2): /* Primitives */
		if( IsPredTok(CurToken) || (CurToken==BackboneTok) )
		{
			if( CurToken == HelixTok || CurToken == SheetTok || CurToken == TurnTok )
			{
				if(pmset != NULL)
				{
					mol_itr = pmset->HostMolecules.begin();
					for( ; mol_itr != pmset->HostMolecules.end(); mol_itr++)
					{
						if( !(*mol_itr)->IsSecStructFound() ) p_mol_editor->DetermineSecStructure((*mol_itr),False);
					}
				}
			}
			if( CurToken == CystineTok)
			{
				if(pmset != NULL)
				{
					if( !pmset->SSBonds_found ) p_mol_editor->FindDisulphideBridges(pmset);
				}
			}

			switch( CurToken )
			{
			case(HelixTok):
				pred = PredHelix;
				break;
			case(SheetTok):
				pred = PredSheet;
				break;
			case(TurnTok):
				pred = PredTurn;
				break;
			case(CystineTok):
				pred = PredCystine;
				break;
			case(BackboneTok):
				pred = PredMainChain;   break;
			case(SelectedTok):
				pred = PropSelect;      break;
			default:
				pred = PredAbsChr(PredTokOrd(CurToken));
			}
			
			tmp1 = new AtomExpr();
			tmp1->type = OpConst|OpLftProp|OpRgtVal;
			tmp1->lft.val = pred;
			FetchToken();
			return( tmp1 );
		}
		else if( IsPropTok(CurToken) )
		{
			int prop_code = 0;
			tmp1 = new AtomExpr();
			tmp1->type = OpLftProp|OpRgtVal;
			switch( CurToken )
			{
				case(TemperatureTok): pred = PropTemp;    break;
				case(RadiusTok):      pred = PropRad;     break;
				case(AtomNoTok):      pred = PropIdent;   break;
				case(ElemNoTok):      pred = PropElemNo;  break;
				case(ChemGrpTok):     pred = PropChemGroup; break;
				case(ResNoTok):       pred = PropResId;   break;
				case(ModelTok):       pred = PropModel;   break;
			}
			tmp1->lft.val = pred;
			prop_code = pred;
			
			FetchToken();
			if( CurToken=='=' )
			{
				tmp1->type |= OpEqual;
				if( FetchToken() == '=' ) FetchToken();
			}
			else if( CurToken=='<' )
			{
				FetchToken();
				if( CurToken=='>' )
				{
					tmp1->type |= OpNotEq;
					FetchToken();
				}
				else if( CurToken=='=' )
				{
					tmp1->type |= OpLessEq;
					FetchToken();
				}
				else tmp1->type |= OpLess;
			}
			else if( CurToken=='>' )
			{
				if( FetchToken() == '=' )
				{
					tmp1->type |= OpMoreEq;
					FetchToken();
				}
				else tmp1->type |= OpMore;
			}
			else if( (CurToken=='!') || (CurToken=='/') )
			{
				if( NextIf('=',"Syntax error in expression") )
				{
					delete tmp1;
					return( (AtomExpr*)NULL );
				}
				else
				{
					tmp1->type |= OpNotEq;
				}
				FetchToken();
			}
			else
			{
				PrintLog("Syntax error in expression\n");
				delete tmp1;
				return( (AtomExpr*)NULL );
			}			
			
			if( CurToken == '-' )
			{
				FetchToken();
				neg = True;
			}
			else
				neg = False;
			
			if( CurToken != NumberTok )
			{
				if( prop_code == PropAtGroup)
				{
					std::string grp_name = TokenIdent;
					AtomGroup* patl = pmset->GetAtomGroupByID(grp_name.c_str());
					if(patl == NULL)
                                        {
						PrintLog(" No atom group with name %s \n",grp_name.c_str());
                                                delete tmp1;
				                return( (AtomExpr*)NULL );
					}
					else
					{
						TokenValue = (int)patl;
					}
				}
				else
				{
				   PrintLog("Integer value expected\n");;
				   delete tmp1;
				   return( (AtomExpr*)NULL );
				}
			}
			
			if( neg )
			{
				tmp1->rgt.val = -(int)TokenValue;
			}
			else
				tmp1->rgt.val = (int)TokenValue;
			FetchToken();
			return( tmp1 );
		}
		else
		{
			switch( CurToken )
			{
			case('('):    FetchToken();
				if( !(tmp1=ParseExpression(0,pmset)) )
					return( (AtomExpr*)NULL );
				
				if( CurToken!=')' )
				{
					PrintLog("Close parenthesis ')' expected\n");
					delete tmp1;
					return( (AtomExpr*)NULL );
				}
				FetchToken();
				return(tmp1);
				
			case('!'): case('~'):
			case(NotTok): FetchToken();
				if( !(tmp1=ParseExpression(2,pmset)) )
					return( (AtomExpr*)NULL );
				
				tmp2 = new AtomExpr();
				tmp2->type = OpNot | OpRgtVal;
				tmp2->lft.ptr = tmp1;
				return( tmp2 );
				
			case('-'):    if( NextIf(NumberTok,"Integer value expected") )
							  return( (AtomExpr*)NULL );
				return( ParseRange(True) );
				
			case(NumberTok):
				return( ParseRange(False) );
				
			case(WithinTok):
				if( NextIf('(',"Open parenthesis '(' expected") )
					return( (AtomExpr*)NULL );


				dtmp =0.0;
				FetchToken();
				if( CurToken==NumberTok )
				{
					dtmp= (double) TokenValue;
				}
				else if( CurToken == FloatTok)
				{
					dtmp = TokenValueFloat;
				}
				else
				{
					PrintLog("Integer or Float value expected\n");;
					return( (AtomExpr*)NULL );
				}
				
				if( dtmp > 50.0 )
				{
					PrintLog("Parameter value too large\n");
					return( (AtomExpr*)NULL );
				}
					
				if( NextIf(',',"Comma separator missing") )
					return( (AtomExpr*)NULL );
				
				FetchToken();
				if( !(tmp1=ParseExpression(0,pmset)) )
					return( (AtomExpr*)NULL );
				
				if( CurToken!=')' )
				{
					PrintLog("Close parenthesis ')' expected\n");
					delete tmp1;
					return( (AtomExpr*)NULL );
				}
				
				FetchToken();
				if( dtmp == 0.0 )
					return( tmp1 );
				
				tmp2 = new AtomExpr();
				tmp2->type = OpWithin;
				dtmp= dtmp;
				tmp2->lft.dval = dtmp*dtmp;
				tmp2->rgt.set = new AtomGroup(tmp1, pmset);
				delete tmp1;
				return( tmp2 );
				
			default:
				if( CurToken==IdentTok )
				{
					tmp1 = AtomExpr::LookUpAtGroupExpr(TokenIdent.c_str(),pmset);
					if( !tmp1 )
						tmp1 = AtomExpr::LookUpElement(TokenIdent.c_str());
					 
					if( tmp1 )
					{
						FetchToken();
						return(tmp1);
					}
				}
				
				cursor_pos = str_start_pos;
				AtomExpr* p_expr = AtomExpr::ParsePrimitiveExpr(pmset,CurLine.c_str(),cursor_pos);
				FetchToken();
				
 				if( p_expr == NULL )
				{
					PrintLog("Syntax error in expression\n");
					return( NULL );
				}
				else
				{
					return( p_expr );
				}
            }
		}
    }
    return( (AtomExpr*)NULL );
}

//#ifdef _WIN32
///* Avoid Optimizer Warning */
//#if ! defined(TWINE)
//    #pragma optimize("g",off)
//#endif
//#endif

