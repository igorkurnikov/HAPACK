#include "parallel.h"
#include "vtype.h"



#include <numer.h>

#include "io_incl.h"
#include "iphys_file.h"


static int ppp;

inline void to_physics(short unsigned int  c1,short unsigned int  c2,
		       short unsigned int  c3,short unsigned int  c4,
		       short unsigned int& p1,short unsigned int& p2,
		       short unsigned int& p3,short unsigned int& p4,
		       double v)
{
  p1 = c1;
  p2 = c3;
  p3 = c2;
  p4 = c4;
  // cout << "XXX" << setw(4) << ppp << setw(4) << p1 << setw(4) << p2 
  // << setw(4) << p3 << setw(4) << p4 << setw(15) << 100*v << endl;

}

inline int swap(int& p1,int& p2)
{
  int tmp;
  if (p1 > p2)
  {
    tmp = p1;
    p1  = p2;
    p2  = tmp;
  }
}

IntegPhysFile::IntegPhysFile(IntegFile& chemfile,int debug) 
  : IntegCoreFile("IPHYS",chemfile.no_orbs(),debug),
    BucketFile<IntegData> (256 + chemfile.no_orbs()*(chemfile.no_orbs()+1),
			   core_file_bucket_size(chemfile.no_orbs()),debug)

{
  const double epsilon_qci = 1e-6;

  if (debug)
    cout << " +++ Entered IntegPhysFile Constructor" << endl;

  int pip,i;
  for(pip = 0; pip < integ_file_twoel_pip_offset; pip++)
  {
    if (!local(pip)) continue;
    const Bucket<IntegData>& b = chemfile.read(pip,0);
    if(b.size() > 0 && debug)
      cout << " +++ Copying PIP: " << setw(4) << pip << " with: " 
	   << setw(8) << b.size() 
	   << " elements on node: " << setw(4) << pm -> id() << endl;
    
    for(i = 0; i < b.size(); i++)
      add(pip,b(i));
  }

  short unsigned int c1,c2,c3,c4; // chemical notation pips
  short unsigned int p1,p2,p_pip;
  
  IntegData ival;
  int mc = 0;
  for(c1 = 0; c1 < no_orbs(); c1++)
    for(c2 = 0; c2 <= c1; c2++)
    {
      pip = integ_file_twoel_pip_offset + c1*(c1+1)/2 + c2;
      if (!local(pip)) continue;
      const Bucket<IntegData>& b = chemfile.read(pip,0);
      if (b.size() > 0 && debug)
	cout << " +++ Sorting 2P: " << setw(3) << c1 << " " << setw(3) << c2 
	     << " PIP: " << setw(4) << pip << " with: " << setw(8) << b.size()
	     << " elements on node: " << setw(4) << pm -> id() << endl;

      ppp = pip - integ_file_twoel_pip_offset;
      
      for(i = 0; i < b.size(); i++)
      {
	c3 = b(i).idx1;
	c4 = b(i).idx2;
	ival.val = b(i).val;
	if (fabs(ival.val) < epsilon_qci){
	  mc++;
	  continue;
	}

	if (c1 == c2 && c1 == c3 && c1 == c4){
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
	else if (c1 == c2 && c3 == c4){
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
	else if (c1 == c2){
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c1,c2,c4,c3,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c4,c3,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
	else if (c3 == c4){
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c2,c1,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c2,c1,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
	else if (c1 == c3 && c2 == c4){
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c2,c1,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c2,c1,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c2,c1,c4,c3,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
	else if (c1 == c4 && c2 == c3){
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c2,c1,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c1,c2,c4,c3,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
	else {
	  to_physics(c1,c2,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c2,c1,c3,c4,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c3,c4,c2,c1,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c1,c2,c4,c3,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c4,c3,c1,c2,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c2,c1,c4,c3,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	  to_physics(c4,c3,c2,c1,p1,p2,ival.idx1,ival.idx2,ival.val);
	  dump(p1,p2,ival);
	}
      }
    }

  cout << "deleted " << mc << " integrals" << endl;
  write_all();  
  
}







