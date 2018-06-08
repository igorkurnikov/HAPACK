/*!  \file habond.h

     Classes to define covalent bond object in HARLEM

    \author Igor Kurnikov
    \date 1997-2002
*/
#ifndef HABOND_H
#define HABOND_H

#include "haconst.h"

class Vec3D;
class HaAtom;
class AtomContainer;

/* Bond Flags */
const int DrawBondFlag =  0x0e;

const int WireFlag     =   0x02;     //!< Depth-cued wireframe         
const int DashFlag     =   0x04;     //!< Dashed Depth-cued wireframe  
const int CylinderFlag =   0x08;     //!< Line/Cylinder representation 

const int HydrBondFlag  =  0x00;     //!< Hydrogen bond [place keeper] 

const int NormBondFlag  =  0x10;
const int DoubBondFlag  =  0x20;
const int TripBondFlag  =  0x40;
const int AromBondFlag  =  0x80;

class HaBond
//!  Class to define Bond object 
{
public:
  HaBond();
  HaBond(const HaBond& ref);
  HaBond(HaAtom* src, HaAtom* dst );
  virtual ~HaBond();

  bool SetParamFrom(const HaBond& ref); 
  bool SetTypeFrom(const HaBond& ref); //!< Set Bond Type from another bond 

  HaBond* assign(HaAtom* src, HaAtom* dst);

  bool operator==(const HaBond & rhs) const;
  bool operator< (const HaBond & rhs) const;

  bool IsSingle()   const { return flag & NormBondFlag; }
  bool IsAromatic() const { return flag & AromBondFlag; }
  bool IsDouble()   const { return flag & DoubBondFlag; }
  bool IsTriple()   const { return flag & TripBondFlag; }
  bool IsVirtual()   const { return type == VIRTUAL_BOND; }

  void SetSingle()   { flag &= 0x0F; flag |= NormBondFlag; } 
  void SetAromatic() { flag &= 0x0F; flag |= AromBondFlag; } 
  void SetDouble()   { flag &= 0x0F; flag |= DoubBondFlag; } 
  void SetTriple()   { flag &= 0x0F; flag |= TripBondFlag; } 
  void SetVirtual()   { type = VIRTUAL_BOND; } 

// Properties Manipulation Utilities:

  void Select();
  void UnSelect();
  int Selected() const;

  void SetNotDraw();  //!< Set bond not to be displayed 
  void DrawWire();   //!< Set to display bond as a wire
  void DrawDashed(); //!< Set to display bond as a dashed line
  void DrawCylinder(double rad); //!< Set to display bond as a cylinder , rad in Bohrs
  int  IsToDraw() const; //!< Check if the bond is set to be displayed
  
// RASMOL BOND structure parameters

  HaAtom* GetFirstAtom()  { return srcatom; } //!< Get first  atom of the bond
  const HaAtom* GetFirstAtom() const  { return srcatom; } //!< Get first  atom of the bond (const version)
  HaAtom* GetSecondAtom() { return dstatom; } //!< Get second atom of the bond
  const HaAtom* GetSecondAtom() const { return dstatom; } //!< Get second atom of the bond

  HaAtom  *srcatom;      //!< Source Atom Ptr       
  HaAtom  *dstatom;      //!< Destination Atom Ptr  

  double radius;         //!< Radius of the bond stick in Bohrs (for images)        
  short irad;            //!< Image Radius          
  short col;             //!< Bond Colour           
  unsigned char  flag;   //!< Database flags        

  enum CovBondType { SINGLE_BOND = 0, DOUBLE_BOND, TRIPLE_BOND, AROMATIC_BOND, VIRTUAL_BOND } type;

protected:

};

typedef std::vector<HaBond*> BondArray;

class HaHBond
//!  Class to define Hydrogen Bond object 
{
public:
  HaHBond();
  HaHBond(HaAtom* new_donor, HaAtom* new_acceptor, HaAtom* h_atom = NULL );
  HaHBond(const HaHBond& ref);
  virtual ~HaHBond();
  
  bool operator==(const HaHBond & rhs) const;
  bool operator< (const HaHBond & rhs) const;

// Flag Manipulation Utilities:

  void  Select();
  void UnSelect();
  int Selected() const;
  
  void SetNotDraw();    //!< Set H-bond not to be displayed 
  void SetDraw();       //!< Set H-bond to be displayed 
  void DrawWire();   //!< Display bond as a wire
  void DrawDashed(); //!< Display bond as a dashed line
  void DrawCylinder(double rad); //!< Display bond as a cylinder , rad in Bohrs
  int IsToDrawCylinder() const; //!< Check of H-bond to be displayed as a cylinder 
  int IsToDrawWire() const; //!< Check of H-bond to be displayed as a wire
  int IsToDrawDashed() const; //!< Check of H-bond to be displayed as a dashed line
  int IsToDraw() const; //!< Check if the bond is set to be displayed

  HaAtom* GetDonorAtom()    { return src; } //!< Get Donor atom of the hydrogen bond
  HaAtom* GetAcceptorAtom() { return dst; } //!< Get accepotr atom of the hydrogen bond
  HaAtom* GetHAtom()        { return p_hatom; } //!< Get accepotr atom of the hydrogen bond
  int GetHCoord(Vec3D& pt);   //!< Get Coordinates of H atom of the donor closest to acceptor 
  
// HBOND structure parameters

  HaAtom*  srcCA;              //!< Source Alpha Carbon      
  HaAtom*  dstCA;              //!< Destination Alpha Carbon 
  HaAtom*  dst;                //!< Acceptor [=CO] Atom Ptr  
  HaAtom*  src;                //!< Donor [=NH] Atom Ptr
  HaAtom*  p_hatom;             //!< Hydrogen atom 

  short energy;                //!< Hydrogen bond energy     
  double radius;               //!< World Radius in Bohrs    
  short irad;                  //!< Image Radius             
  char offset;                 //!< Signed difference between residue numbers of the atoms forming hydrogen bond         
  short col;                   //!< Hydrogen bond colour     
  
protected:

  unsigned char flag;          //!< Database flags

};

#endif /* !HABOND_H */
