#include <numer.h>
#include "vtype.h"
#include "integ_file.h"

int integ_file_max(int a,int b) { return (a < b) ? b : a; }
int core_file_bucket_size(int no_orbs) 
{
  return integ_file_max(10,no_orbs*(no_orbs+1)/2/100);
}

/*TEX
\subsection{Constructors}
*/
IntegFile::IntegFile(int no_orbs) :
  BucketFile<IntegData> (256 + no_orbs*(no_orbs+1),
			 integ_file_max(20,no_orbs*(no_orbs+1)/2/100),0)
{
  orbital_map.reset(no_orbs);
  reset_orbital_map();
}

void IntegFile::reset_orbital_map()
{
  for(int i = 0; i < orbital_map.size(); i++)
    orbital_map[i] = i;
}

void IntegFile::sp_hamiltonian(IntegMat& m)
{
  IntegMat tmp;
  read_mat(m,IntegKinetic);
  read_mat(tmp,IntegNuclear);
  m -= tmp;
  read_mat(tmp,IntegECP);
  m += tmp;
  read_mat(tmp,IntegCore);
  m += tmp;
}

IntegDiskFile::IntegDiskFile(const char* fname,int no_orbs,int debug) :
  BucketFile<IntegData> (256 + no_orbs*(no_orbs+1),
			 integ_file_max(10,no_orbs*(no_orbs+1)/2/100),
			 debug),
  IntegFile(no_orbs),
  FlexFile<IntegData> (fname,256 + no_orbs*(no_orbs+1),
		       integ_file_max(20,20 + no_orbs*(no_orbs+1)/2/100),debug)

{
  data.idx1 = 0;
  data.idx2 = 0;
  data.val  = no_orbs;
  add(0,data);
}

IntegCoreFile::IntegCoreFile(const char* fname,int no_orbs,int debug) :
  IntegFile(no_orbs),
  CoreFile<IntegData> (256 + no_orbs*(no_orbs+1),
		       core_file_bucket_size(no_orbs),debug),
  BucketFile<IntegData> (256 + no_orbs*(no_orbs+1),
			 core_file_bucket_size(no_orbs),debug)
{
  data.idx1 = 0;
  data.idx2 = 0;
  data.val  = no_orbs;
  add(0,data);
}
 

IntegDiskFile::IntegDiskFile(const char* fname) : 
  BucketFile<IntegData> (256,12,0),
  IntegFile(0), FlexFile<IntegData> (fname)
{
  const Bucket<IntegData>& reader = FlexFile<IntegData>::read(0,0);
  int no_orbs = (int) (reader(0).val + 0.5);
  orbital_map.reset(no_orbs);
  reset_orbital_map();
}



IntegDiskFile::~IntegDiskFile() {}
IntegCoreFile::~IntegCoreFile() {}
IntegFile::~IntegFile() {}
 
void IntegFile::write_mat(IntegMat& mat,int target_key)
{
  if (!pm->master())
    return;
  read(target_key,1); // clean up first
  
  for(int i = 0; i <  mat.size1(); i++)
  for(int j = 0; j <= i; j++)
  {
    if (fabs(mat(i,j)) < epsilon)
	continue;
    data.idx1 = orbital_map(i);
    data.idx2 = orbital_map(j);
    data.val = mat(i,j);
    add(target_key,data);
  }
}

int IntegFile::read_mat(IntegMat& mat,int target_key,int clear_code)
{
  if (local(target_key))
  {
    int no_orbs = orbital_map.size();
    const Bucket<IntegData>& reader = read(target_key,clear_code);
    mat.reset(no_orbs,no_orbs); 
    mat.set(0);

    
    int o1,o2;
    for(int i = 0; i < reader.size(); i++)
    {
      o1 = reader(i).idx1;
      o2 = reader(i).idx2;
      mat(o1,o2) = reader(i).val;
      mat(o2,o1) = reader(i).val;
    }
  }
  //cout << " Node: " << pm -> id() << " broadcasting mat: " << target_key 
  //     << " with addr: " << addr(target_key) << endl;
  pm -> broadcast(mat,addr(target_key)); // get the data
  double sum = 0;
  for(int i = 0; i < mat.size1(); i++)
    for(int j = 0; j < mat.size2(); j++)
      sum += fabs(mat(i,j));
  if (sum > 0)
    return 1;
  return 0;
  
}

void ifile_gnu_dummy()
{
  ARRAY<IntegData> ia;
  ARRAY<IntegData> ib;
  ib = ia;
}


/*TEX
\subsection{Output Operators}
*/

ostream& operator<<(ostream& os,IntegData data)
{
  char buf[100];
  sprintf(buf," %10.6f",data.val);
  return os << setw(4) << data.idx1 << setw(4) << data.idx2 << buf;
}

  
