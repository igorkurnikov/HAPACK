class DataVault;

class Element
{
  int cc;  
  static int iflag;
  static ARRAY<String> ename;
  static ARRAY<int>    coresize;
friend  istream& operator>>(istream& is,Element& el);

public:
       Element() { cc = 0; } 
       Element(int chg);
       Element(String& nme);
      ~Element() {}
String name()   const { return ename(cc);    }     
int    charge() const { return cc;           }   
int    core()   const { return coresize(cc); } 
static void init();  
static void clean();  
};

istream& operator>>(istream& is,Element& el);
inline ostream& operator<<(ostream& os,const Element el)
{ return os << el.name(); }


class Molecule  
{
           SStack<Location,max_cntr>  ctr;
           SStack<Element ,max_cntr>  elem;
           SStack<String  ,max_cntr>  bname;
public:
           Molecule() {}  
          ~Molecule() {}  
void       add(const Element& el,const Location& loc,const String& bn)
           { ctr.add(loc); elem.add(el); bname.add(bn); }
int        size() const  { return ctr.size(); } 
Location   center (const int i) const { return ctr (i); } 
Element    element(const int i) const { return elem(i); }
const String& basis(const int i) const { return bname(i); }
void       write(DataVault& vault);
void       read (DataVault& vault);
};

ostream& operator<<(ostream& os,const Molecule& molec);
ostream& operator>>(ostream& os,const Molecule& molec);

