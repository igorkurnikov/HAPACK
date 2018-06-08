#include <numer.h>
#include <errorip.h>

class BiSpline
{
private:
  int dimx;
  int dimy;

  D_Array<D_Vector> data;
  D_Array<D_Vector> y2a;
  D_Vector x1a;
  D_Vector x2a;

  D_Vector ytmp;
  D_Vector yytmp;
  
void   spline(D_Vector& x,D_Vector& y,int n,double yp1,double ypn,D_Vector& y2);
void   splie2();
double splint(D_Vector& xa,D_Vector& ya,D_Vector& y2a,double x);

void   init(); // after Data, xa1 and xa2 have been initialized build tables

public:

       BiSpline(char* fname);
       BiSpline(int testfun);

double BiSpline::evaluate(double x1,double x2);
double minx() const { return x1a(   1); }
double maxx() const { return x1a(dimx); }
double miny() const { return x2a(   1); }
double maxy() const { return x2a(dimy); }
  
};

BiSpline::BiSpline(char* fname)
{
  ifstream infile(fname,ios::in);
  if (!infile)
    error(" File not found. ");

  // read x and y data points
  
  int i,j;
  infile >> dimx;
  if (!infile)
    error(" Read Error on DIMX.");
  cout << " # of X Grind Points: " << dimx << endl;
  x1a.reset(dimx,1);
  
  for(i = 1; i <= dimx; i++)
  {
    infile >> x1a[i];
    if (!infile)
      error(" Read Error on XVAL.", i);
  }
  
  infile >> dimy;
  if (!infile)
    error(" Read Error on DIMY.");
  cout << " # of Y Grind Points: " << dimy << endl;
  x2a.reset(dimy,1);
  
  for(i = 1; i <= dimy; i++)
  {
    infile >> x2a[i];
    if (!infile)
      error(" Read Error on YVAL.", i);
  }
  
  // read data

  double d1,d2,val;
  data.reset(dimx,1);
  double dx = (x1a[dimx] - x1a[1])/dimx;
  double dy = (x2a[dimy] - x2a[1])/dimy;
  
  for(i = 1; i <= dimx; i++)
  {
    data[i].reset(dimy,1);
    for(j = 1; j <= dimy; j++)
    {
      infile >> d1 >> d2 >> val;
      int ii = (int) ((d1 - x1a[1])/dx) + 1;
      int jj = (int) ((d2 - x2a[1])/dy) + 1;
      cout << d1 << " " << d2 << " " << ii << " " << jj;
      data[ii][jj] = val;
      cout << "  " << data[ii][jj] << endl;
      if (!infile)
	error(" Read Error on Data.");
    }
  }
  
  init();

}

BiSpline::BiSpline(int test_fun)
{
  int dimx = 10;
  int dimy = 20;
  
  data.reset(dimx,1);
  x1a.reset (dimx,1);
  x2a.reset (dimy,1);
  
  for(int ix = 0; ix < dimx; ix++)
  {
    double x = (double) ix / (dimx - 1);    
    x1a [1 + ix] = x;
    data[1 + ix].reset(dimy,1);
    
    for(int iy = 0; iy < dimy; iy++)
    {
      double y = (double) iy / (dimy - 1);    
      x2a [iy + 1] = y;      
      data[ix + 1][iy + 1] = (2 + x)*(1.7 - y);
    }
  }
  
  init();

  for(int i = 0; i < 100; i++)
  {
    double x = drand48();
    double y = drand48();
    double f  = evaluate(x,y);
    double ff = (2 + x)*(1.7-y);
    if (fabs(f - ff) > 1E-3)
      cout << " Error for: " << x << " " << y << " " << f << " " << ff << endl;
  }
}


void BiSpline::init()
{
  dimx = data.size();
  dimy = data[1].size();
  
  ytmp .reset(dimx,1);
  yytmp.reset(dimx,1);
  y2a.reset(dimx,1);

  for (int j=1;j <= dimx;j++)
  {
    y2a[j].reset(dimy,1);
    spline(x2a,data[j],dimy,1.0e30,1.0e30,y2a[j]);
  }
  
}

void BiSpline::spline(D_Vector& x,D_Vector& y,int n,double yp1,double ypn,D_Vector& y2)
{
  int i,k;
  double p,qn,sig,un;

  D_Vector u(n-1,1);
  
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else 
  {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) 
  {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

double BiSpline::evaluate(double x1,double x2)
{
  for (int j=1;j<=dimx;j++)
    yytmp[j] = splint(x2a,data[j],y2a[j],x2);
  spline(x1a,yytmp,dimx,1.0e30,1.0e30,ytmp);
  double y = splint(x1a,yytmp,ytmp,x1);
  return y;
}

double BiSpline::splint(D_Vector& xa,D_Vector& ya,D_Vector& y2a,double x)
{
  int klo,khi,k;
  double h,b,a;
  int n = xa.size();
  
  klo=1;
  khi=n;
  while (khi-klo > 1) 
  {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) 
    error("Bad xa input to routine splint");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  double y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return y;
}


class Triangulate : public BiSpline
{
  
public:
  Triangulate(char* infile,int nx,int ny);
};

class Mesh : public BiSpline
{
  
public:
  Mesh(char* infile,int nx,int ny);
};

Triangulate::Triangulate(char* infile,int nx,int ny) : BiSpline(infile)
{
  double dx = (maxx() - minx())/nx;
  double dy = (maxy() - miny())/ny;
  
  ofstream outfile("data.pov",ios::out);
  
  // outfile << "$ DATA=CURVE3D\n% dfilltype=1\n%linetype=0\n% pointid=False\n% hiddenline=True" << endl;

  int ix,iy;

  double x00 = minx();
  double y00 = miny();
  double f00 = evaluate(x00,y00);
  double maxf = f00;
  double minf = f00;
  
  for(ix = 0; ix <= nx; ix++) 
    for(iy = 0; iy <= ny; iy++) 
    {
      double x00 = minx() + ix*dx;
      double y00 = miny() + iy*dy;
      double f00 = evaluate(x00,y00);
      if (f00 > maxf) maxf = f00;
      if (f00 < minf) minf = f00;
    }

  cout << maxf << " " << minf << endl;
  
  for(ix = 0; ix < nx; ix++) 
    for(iy = 0; iy < ny; iy++) 
    {
      double x00 = minx() + ix*dx;
      double y00 = miny() + iy*dy;
      double f00 = evaluate(x00,y00);
      
      double x10 = x00 + dx;
      double y10 = y00;
      double f10 = evaluate(x10,y10);
      
      double x01 = x00;
      double y01 = y00 + dy;
      double f01 = evaluate(x01,y01);

      double x11 = x00 + dx;
      double y11 = y00 + dy;
      double f11 = evaluate(x11,y11);
      
      int c = (int) (5 * ((f00 + f01 + f10)/3 - minf) / (maxf-minf));
      // outfile << "\n% fillcolor=" << c << endl;
      outfile << "triangle { < " << x00 << "," << f00 << "," << y00 << "> ," << endl; 
      outfile << "           < " << x01 << "," << f01 << "," << y01 << "> ," << endl; 
      outfile << "           < " << x10 << "," << f10 << "," << y10 << ">  texture { TRed } }" << endl;

      outfile << "triangle { < " << x11 << "," << f11 << "," << y11 << "> ," << endl; 
      outfile << "           < " << x01 << "," << f01 << "," << y01 << "> ," << endl; 
      outfile << "           < " << x10 << "," << f10 << "," << y10 << ">  texture { TRed }}" << endl;

      /*
	c = (int) (5 * ((f11 + f01 + f10)/3 - minf) / (maxf-minf));
      outfile << "\n% fillcolor=" << c << endl;
      outfile << x01 << " " << y01 << " " << f01 << endl; 
      outfile << x10 << " " << y10 << " " << f10 << endl; 
      outfile << x11 << " " << y11 << " " << f11 << endl; 
      */
  }
  
}

Mesh::Mesh(char* infile,int nx,int ny) : BiSpline(infile)
{
  double dx = (maxx() - minx())/nx;
  double dy = (maxy() - miny())/ny;
  
  ofstream outfile("mesh",ios::out);
  ofstream outx("xx",ios::out);
  ofstream outy("yy",ios::out);
  
  // outfile << "$ DATA=CURVE3D\n% dfilltype=1\n%linetype=0\n% pointid=False\n% hiddenline=True" << endl;

  int ix,iy;

  double x00 = minx();
  double y00 = miny();
  double f00 = evaluate(x00,y00);
  double maxf = f00;
  double minf = f00;

  for(ix = 0; ix <= nx; ix++) 
    outx << minx() + ix*dx << " ";
  outx << endl;

  for(iy = 0; iy <= ny; iy++) 
    outy << miny() + iy*dy << " ";
  outy << endl;
  
  for(ix = 0; ix <= nx; ix++) 
  {
    for(iy = 0; iy <= ny; iy++) 
    {
      double x00 = minx() + ix*dx;
      double y00 = miny() + iy*dy;
      double f00 = evaluate(x00,y00);
      outfile << f00 << " ";
    }
    outfile << endl;
  }
  
  return;
  
  for(ix = 0; ix < nx; ix++) 
    for(iy = 0; iy < ny; iy++) 
    {
      double x00 = minx() + ix*dx;
      double y00 = miny() + iy*dy;
      double f00 = evaluate(x00,y00);
      
      double x10 = x00 + dx;
      double y10 = y00;
      double f10 = evaluate(x10,y10);
      
      double x01 = x00;
      double y01 = y00 + dy;
      double f01 = evaluate(x01,y01);

      double x11 = x00 + dx;
      double y11 = y00 + dy;
      double f11 = evaluate(x11,y11);
      
      int c = (int) (5 * ((f00 + f01 + f10)/3 - minf) / (maxf-minf));
      // outfile << "\n% fillcolor=" << c << endl;
      outfile << "triangle { < " << x00 << "," << f00 << "," << y00 << "> ," << endl; 
      outfile << "           < " << x01 << "," << f01 << "," << y01 << "> ," << endl; 
      outfile << "           < " << x10 << "," << f10 << "," << y10 << ">  texture { TRed } }" << endl;

      outfile << "triangle { < " << x11 << "," << f11 << "," << y11 << "> ," << endl; 
      outfile << "           < " << x01 << "," << f01 << "," << y01 << "> ," << endl; 
      outfile << "           < " << x10 << "," << f10 << "," << y10 << ">  texture { TRed }}" << endl;

      /*
	c = (int) (5 * ((f11 + f01 + f10)/3 - minf) / (maxf-minf));
      outfile << "\n% fillcolor=" << c << endl;
      outfile << x01 << " " << y01 << " " << f01 << endl; 
      outfile << x10 << " " << y10 << " " << f10 << endl; 
      outfile << x11 << " " << y11 << " " << f11 << endl; 
      */
  }
  
}


main(int argc,char** argv)
{
  int dimx = 10;
  int dimy = 20;
  
  Mesh t(argv[1],atoi(argv[2]),atoi(argv[3]));
  
}
