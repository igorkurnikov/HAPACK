#include <iostream.h>
#include <stdio.h>
#include <iomanip.h>
#include <string.h>
#include <xtime.h>

enum TNames  { tmain, tsub1,tsub2,tsub3,tsub11,tsub12,tsub13 };
char* tstring[] = { "TMAIN" , "TSUB1" , "TSUB2" , "TSUB3" , "TSUB11" , 
                   "TSUB12", "TSUB13" };  
  
int  count[100];
int  units[100];

double use(int c)
{
  double x = 1;
  
  for(int i = 0; i <     c; i++)
  for(int j = 0; j < 100000; j++)
    x += 2.3*j;
  return x;
}

double sub11(int k)
{
  timer -> start(tsub11);
  count[tsub11]++;
  units[tsub11] += k;
  units[ tsub1] += k;
  units[tmain]  += k;
  double x = use(k);
  timer -> stop(tsub11);
  return x; 
}

double sub12(int k)
{
  timer -> start(tsub12);
  count[tsub12]++;
  units[tsub12] += k;
  units[ tsub1] += k;
  units[tmain]  += k;
  double x = use(k);
  timer -> stop(tsub12);
  return x;
}


double sub1(int k)
{
  timer -> start(tsub1);
  count[tsub1]++;
  units[tsub1] += 20*k;
  units[tmain] += 20*k;
  double x = use(10*k);
  x += sub11(10*k);
  x += use(10*k);
  for(int i = 0; i< 5*k; i++)
    x += sub12(i+3);

  timer -> stop(tsub1);
  return x;
}

main(int argc,char** argv)
{
  cout << "size: " << sizeof(clock_t) << endl;
  
  int cc = 1;
  
  if (argc > 1) 
    sscanf(argv[1],"%i",&cc);
  cout << "Count: " << cc << endl;
  
  timer = new Timer(20);

  int i;
  for(i = 0; i <= tsub13; i++)
  {
    timer -> set_name(i,tstring[i]);
    count[i] = 0;
    units[i] = 0;
  }
  
  timer -> start(tmain);
  count[tmain] += 1;
  units[tmain] += 30;
  double x = use(30);
  x += sub1(cc);
  
  timer -> stop(tmain);
  timer -> print();

  cout << "EXPLICIT COUNT: " << x <<  endl;
  
  double utime = ((double) clock()) / units[tmain] / CLOCKS_PER_SEC*1000.0;
  cout << "MSEC/UNIT: " << utime << endl;
  int    ops = 0;
  
  for(i = 0; i <= tsub13; i++)
  {
    cout << setw(20) << tstring[i] << setw(10) << count[i] << "  " 
      << setw(10) << units[i] << "  " << setw(10) << (int) (units[i]*utime) 
      << endl;
    ops += units[i];
  }

  cout << "Flop Rate: " << ops * 2E5 / clock()  << endl;
  
}
