/*TEX
\subsection{Utility for Symmetry Adaptation} 

*/

#include "io_incl.h"
#include <numer.h>
#include "typesip.h"
#include <tsstack.h>

void model(SStack<int>& wave)
{
  Momentum m;
  int count = 0;
  for(int i=0; i<wave.size(); i++)
    count += m.number(wave[i]);

  cout << "A total of " << count << " orbitals " << endl;
  
  Matrix t(2*count,2*count);
  t.set(0);
  double root = 1.0/sqrt(2.0);
  
  int start = 0;
  for(i=0; i<wave.size(); i++)
  {
    for(Momentum m(wave[i]); m(); ++m)
    {
      int sgn = (m.x() % 2 == 0) ? 1 : -1;
      int idx = start + m.local_index();
      // sgn = 1;
      t(idx      ,      idx)  =      root;
      t(idx      ,count+idx)  = -sgn*root;
      t(count+idx,      idx)  =      root;
      t(count+idx,count+idx)  = +sgn*root;
    }    
    start += m.number(wave[i]);
  }
  
  t.print(cout);
  t.write("TRANS");
}

void smodel(SStack<int>& wave)
{
  int even[4][10];
  
  even[0][0] =  1;

  even[1][0] = -1;
  even[1][1] =  1;
  even[1][2] =  1;

  even[2][0] = -1;
  even[2][1] = -1;
  even[2][2] =  1;
  even[2][3] =  1;
  even[2][4] =  1;

  even[3][0] = -1;
  even[3][1] = -1;
  even[3][2] =  1;
  even[3][3] =  1;
  even[3][4] = -1;
  even[3][5] =  1;
  even[3][6] =  1;

  Momentum m;
  int count = 0;
  for(int i=0; i<wave.size(); i++)
    count += 2*wave[i]+1;

  cout << "A total of " << count << " orbitals " << endl;
  
  Matrix t(2*count,2*count);
  t.set(0);
  double root = 1.0/sqrt(2.0);
  
  int start = 0;
  for(i=0; i<wave.size(); i++)
  {
    for(int m=0; m < 2*wave[i]+1; ++m)
    {
      int sgn = even[wave[i]][m];
      int idx = start + m;
      // sgn = 1;
      t(idx      ,      idx)  =      root;
      t(idx      ,count+idx)  = -sgn*root;
      t(count+idx,      idx)  =      root;
      t(count+idx,count+idx)  = +sgn*root;
    }    
    start += 2*wave[i]+1;
  }
  t.print(cout);
  t.write("TRANS");
}

main()
{
  int sphereflag;
  
  SStack<int> wave;
  
  cin >> sphereflag;
  
  if (sphereflag) cout <<  "SPHERICAL FUNCTIONS" << endl;
  
  while(1)
  {
    int m;
    cin >> m;
    if (m >= 0) 	wave.add(m);
    else break;
  }
  
  cout << " Momenta: " << wave << endl;
  
  if (sphereflag)
    smodel(wave);
  else
    model(wave);
  

}
