/***************************************************************************
 *   Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 2006                                      *
 ***************************************************************************/

#if defined(_MSC_VER)
#include <process.h>
#endif

#include <mpi.h>

#include "command.h"
#include "stdlib.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include "apbsmod.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "etcoupl.h"
#include "hasurface.h"
#include "hamolview.h"
#include "harlemapp.h"

#include "haconsts.h"

double APBSMod::d = 0.0;

APBSMod::APBSMod(MolSet* new_phost_mset) :
HaCompMod(COMP_MOD_ELECTROST,new_phost_mset)
{
  epsi=2.0f;
  epsout=81.0f;
  nx=ny=nz=65;
  GridScale=1.0f;
  rionst=0.0;
  exrad=2.0;
  radprb= (float) 1.4;
  boundary=0;
  iper[0]=0;
  iper[1]=0;
  iper[2]=0;
  nlit=100;
  nnit=0;
  chgm[0]='s';chgm[1]='p';chgm[2]='l';chgm[3]='0';chgm[4]='\0';
  bcfl[0]='z';bcfl[0]='e';bcfl[0]='r';bcfl[0]='o';bcfl[4]='\0';
  PBEType[0]='l';PBEType[0]='p';PBEType[0]='b';PBEType[0]='e';PBEType[4]='\0';
  nlev=4;
}

APBSMod::~APBSMod()
{
}
int APBSMod::SavePQRFile(const char* filename)
{
  if(phost_mset == NULL)return EXIT_FAILURE;
  int atonnum=0,resnum=0;
  double x, y, z,q,r;
  
  HaChain  *chain;
  HaAtom  *aptr;
  int count;

  FILE* DataFile = fopen( filename, "w" );
  if( !DataFile )
  {
    PrintLog("\n");
    return( False );
  }
  
  AtomIteratorMolSet aitr(phost_mset);                                  // FIX IGOR
  for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())      // FIX IGOR
  {
    HaResidue* group = aptr->GetHostRes();                             // FIX IGOR
    if(phost_mset->save_opt_default.save_selected && !aptr->Selected() )continue;
        
    x = aptr->GetX_Ang();
    y = aptr->GetY_Ang();
    z = aptr->GetZ_Ang();
    r = aptr->radius;
    q = aptr->GetCharge();
    
    if( aptr->flag&HeteroFlag )
      fputs("HETATM",DataFile);
    else
      fputs("ATOM  ",DataFile);
    
	std::string atname(aptr->GetName());
    if(atname.size() < 4)
      atname.insert(0," ");

    int k;
    for(k=0; k < 3; k++)
    {
      if(atname.size() < 4)
        atname+= " ";
    }
						
    if(atname.size() > 4)
      atname = atname.substr(0,3);
    
	std::string res_name(group->GetName());

    for(k=0; k < 3; k++)
    {
      if(res_name.size() < 3)
        res_name.insert(0," ");
    }
						
    if(res_name.size() > 3)
      res_name = res_name.substr(0,2);
    
    fprintf( DataFile, "%5d %.4s %.3s %4d    ",atonnum++, atname.c_str(), res_name.c_str(),group->serno );
    fprintf( DataFile, "%8.3f%8.3f%8.3f%8.3f%8.3f\n",x,y,z,q,r);
  }
  fclose( DataFile );
  return EXIT_SUCCESS;
}
int APBSMod::Run()
{
  if(phost_mset == NULL)return EXIT_FAILURE;
  //prepare input files
  //xyzqr
  char PQRFileName[]="apbs.pqr";
  SavePQRFile(PQRFileName);
  //input
  char FileName[]="apbs.in";
  FILE *OUT=fopen(FileName,"w");
  fprintf(OUT,"\n");
  fprintf(OUT,"# READ IN MOLECULES\n");
  fprintf(OUT,"read \n");
  fprintf(OUT,"    mol pqr %s       # Read molecule 1\n",PQRFileName);
  fprintf(OUT,"end\n");
  
  fprintf(OUT,"# CALCULATE POTENTIAL FOR SOLVATED STATE\n");
  fprintf(OUT,"elec name solvated\n");
  fprintf(OUT,"    mg-manual\n");
  fprintf(OUT,"    dime %d %d %d\n",nx,ny,nz);
  fprintf(OUT,"    nlev %d\n",nlev);
  float grd=1.0f/GridScale;
  fprintf(OUT,"    grid %f %f %f\n",grd,grd,grd);
  fprintf(OUT,"    gcent mol 1\n");
  fprintf(OUT,"    mol 1 \n");
  fprintf(OUT,"    %s\n",PBEType);
  fprintf(OUT,"    bcfl %s\n",bcfl);
  fprintf(OUT,"    ion 1.0 %f %f\n",rionst,exrad);
  fprintf(OUT,"    ion -1.0 %f %f\n",rionst,exrad);
  fprintf(OUT,"    pdie %f\n",epsi);
  fprintf(OUT,"    sdie %f\n",epsout);
  fprintf(OUT,"    chgm %s\n",chgm);
  fprintf(OUT,"    srfm mol\n");
  fprintf(OUT,"    srad %f\n",radprb);
  fprintf(OUT,"    swin 0.3\n");
  fprintf(OUT,"    sdens 10.0\n");
  fprintf(OUT,"    temp 298.15\n");
  fprintf(OUT,"    gamma 0.105\n");
  fprintf(OUT,"    calcenergy total\n");
  fprintf(OUT,"    calcforce no\n");
  fprintf(OUT,"end\n");
  
  fprintf(OUT,"# COMBINE TO FORM SOLVATION ENERGY\n");
  fprintf(OUT,"print energy solvated end\n");

  fprintf(OUT,"# SO LONG\n");
  fprintf(OUT,"quit\n");
  
  fclose(OUT);
  
  //Run APBS
  std::string cmdrun="apbs ";
  cmdrun+=FileName;
  cmdrun+=">apbs.out";
  system(cmdrun.c_str());
  //ReadResults
  char line[1024];
  OUT=fopen("apbs.out","r");
  while(!feof(OUT))
  {
    fgets(line,1024,OUT);
    // printf(line);
	std::string HaLine(line);
    if(HaLine.find("Global net energy")!= std::string::npos)
    {
      sscanf(line,"  Global net energy = %lf",&E);
			E=E*HaConsts::kJ_mol_to_kT;
      printf("Cool I find E=%lf kT\n",E);
    }
  }
  fclose(OUT);
  return EXIT_SUCCESS;
}











