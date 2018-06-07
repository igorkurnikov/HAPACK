#include <stringc.h>
#include <tsstack.h>
#include <errorip.h>
#include "typesip.h"
#include "vault.h"
#include "element.h"

static char *element_names[] = { "X","H","He",
   "Li","Be","B","C","N","O","F","Ne",
   "Na","Mg","Al","Si","P","S","Cl","Ar",
   "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As",
   "Se","Br","Kr" }; 
          
ARRAY<String> Element::ename; 
ARRAY<int>    Element::coresize;
int           Element::iflag = 0;  // initcates intialization of arrays;

void Element::clean()
{
  ename.reset(0);
  coresize.reset(0);
  iflag = 0;
}

void Element::init()
{
  int noel = 36;
  ename   .reset(noel+1);
  coresize.reset(noel+1);
  int i;
  for(i =0; i <= noel; i++)
    ename[i] = element_names[i];
  
  coresize[0]  =   0;  
  coresize[1]  =   0;  
  for(i= 2; i <= 10; i++)
    coresize[i]  =     2; 
  for(i=11; i <= 18; i++)
    coresize[i]  =    10;  
  for(i=19; i <= 36; i++)
    coresize[i]  =    18;
  
/*
  for(i = 0 ; i <= noel; i++)
    cout << "Element: " << i << "  " << ename[i] << " " << coresize[i] << endl;
*/
  iflag = 1;
}

Element::Element(int chg)
{
  if (!iflag)
    init();
  
  if (chg >=0 && chg < ename.size())
    cc = chg;
  else
  {
    cout << " +++ Fatal Error in Element::Element -- charge was: " << chg 
         << endl;
    error("Fatal Error.");
  }
  
}

Element::Element(String& nme)
{
  if (!iflag)
    init();
  cc = -1;
  for(int i = 0; i < ename.size(); i++)
    if (nme == ename[i])
      cc = i;
  if (cc == -1)
  {
    cout << " +++ Fatal Error in Element::Element -- name was: " << nme 
         << endl;
    error("Fatal Error.");
  }
}

istream& operator>>(istream& is,Element& el)
{
  String s;
  is >> s;
  if (!is)
    error("IO Error for Element. ");
  Element tmp(s);
  el = tmp;
  return is;
}

ostream& operator<<(ostream& os,const Molecule& molec)
{
  os << molec.size() << " ";
  for(int i = 0; i < molec.size(); i++)
    os << molec.element(i) << " " << molec.center(i) << " ";
  return os;
}

void Molecule::write(DataVault& vault)
{
  NUMARRAY<int> charge(ctr.size());
  Vec           xcoord(ctr.size());
  Vec           ycoord(ctr.size());
  Vec           zcoord(ctr.size());
  for(int i = 0; i < ctr.size(); i++)
  {
    charge[i] = element(i).charge();
    xcoord[i] = ctr(i).x();
    ycoord[i] = ctr(i).y();
    zcoord[i] = ctr(i).z();
  }
  vault.insert("GEOM.CHARGE",charge);
  vault.insert("GEOM.XCOORD",xcoord);
  vault.insert("GEOM.YCOORD",ycoord);
  vault.insert("GEOM.ZCOORD",zcoord);
}



