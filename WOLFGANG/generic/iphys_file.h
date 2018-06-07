#include "integ_file.h"

class IntegPhysFile : public IntegCoreFile
{
  void dump(int p1,int p2,IntegData ival)
    {
      if (p1 >= p2)
      {
	int p_pip    = p1*(p1+1)/2 + p2 + integ_file_twoel_pip_offset;
	add(p_pip,ival);
	if (p1 == ival.idx1 && p2 == ival.idx2)
	   add(IntegDiagDir,ival);
	 if (p1 == ival.idx2 && p2 == ival.idx1)
	   add(IntegDiagExc,ival);
      }
    }
  
public:  
  IntegPhysFile(IntegFile& chemfile,int debug = 0);
  void get_pip(int o1,int o2,ARRAY<IntegData>& iarray)
    {
      int p_pip = (o1 > o2) ? o1*(o1+1)/2 + o2 : o2*(o2+1)/2 + o1;
      p_pip += integ_file_twoel_pip_offset;
      
      const Bucket<IntegData>& b = read(p_pip);
      
      int n = b.size();
      pm -> broadcast(n,addr(p_pip));
      iarray.reset(n);
      if (local(p_pip))
      {
	for(int i = 0; i < n; i++)
	  iarray[i] = b(i);
      }
      int sz = n * sizeof(IntegData);
      pm -> broadcast(iarray.get_base(),sz,addr(p_pip));
    }
  int get_pip_local(int o1,int o2,ARRAY<IntegData>& iarray)
    {
      int p_pip = (o1 > o2) ? o1*(o1+1)/2 + o2 : o2*(o2+1)/2 + o1;
      p_pip += integ_file_twoel_pip_offset;

      if (!local(p_pip)) return 0;
      
      const Bucket<IntegData>& b = read(p_pip);
      int n = b.size();
      iarray.reset(n);
      for(int i = 0; i < n; i++)
	iarray[i] = b(i);
      return 1;
    }
  int no_integrals()
    {
      int c=0;
      for (int i=0;i<number_of_keys();i++)
	if (local(i)){
	  const Bucket<IntegData>& b = read(i);
	  c += b.size();
	}
      return c;
    }
};
  





