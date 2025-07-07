#include "io_incl.h"
#include <errorip.h>
#include "vtype.h"
#include "vault.h"

DataVault::DataVault(String fname,DataVaultState code) : filename(fname)
{
  if (code == vault_old)
    read();
}

DataVault::~DataVault() { write(); } 

void DataVault::write()
{
  char * fn = filename();
  ofstream os(fn,ios::out);
  delete fn;
  if (!os)
  {
    cout << " +++ DataVault could not open: " << filename << endl;
    error("DataVault could not open file. ");
  }
  write(os);
}

void DataVault::read()
{
  char * fn = filename();
  ifstream is(fn,ios::in);
  delete fn;
  if (!is)
  {
    cout << " +++ DataVault could not open: " << filename << endl;
    error("DataVault could not open file. ");
  }
  read(is);
}
  
void DataVault::read(istream& is)
{
  int no_item;
  is >> no_item;
  // cout << " +++ Opening DataVault with " << no_item << " items." << endl;
  nname.resize(no_item+100);
  ncont.resize(no_item+100);
  
  String nn,cc;
  DataVaultType type;
  for(int i = 0; i < no_item; i++)
  {
    is >> type;
    nn.read(is);
    cc.read(is);
    if (!is)
    {
      cout << " +++ DataVault::Severe Error in Reading item: " << i << endl;
      error("DataVault::Severe Error. ");
    }
    insert(nn,cc,type,vault_severe);
  }
}

void DataVault::write(ostream& os)
{
  int no_item = nname.size();
  cout << " +++ Writing DataVault with " << no_item << " items." << endl;
  os << no_item << endl;
  String nn,cc;
  
  for(int i = 0; i < no_item; i++)
  {
    os << ntype[i] << " ";
    nname[i].write(os);
    ncont[i].write(os);
    os << endl;
    if (!os)
    {
      cout << " +++ DataVault::Severe Error in Writing item: " << i << endl;
      error("DataVault::Severe Error. ");
    }
  }
}

int  DataVault::check(const String& name)
{
  return (find(name) >= 0);
}

int  DataVault::check(const String& name,DataVaultType type)
{
  int pos = find(name);
  return (pos >=0) && ntype[pos] == type;  
}

int  DataVault::find(const String& name,DataVaultType type,
		     DataVaultCode code)
{
  for(int i = 0; i < nname.size(); i++)
    if (name == nname[i])
    {
      if (type != v_string && ntype[i] != type)
      {
	cout << " +++ DataVault Error:  Wrong type for item: " 
	     << name << " expected: " << type << " found: " << ntype[i]
	     << endl;
        error("Fatal Error in DataVault. ");
      }
      return i;
    }

  switch(code)
  {
  case vault_ignore:
    break;
  case vault_warning:
    cout << " +++ DataVault Warning: Did not find item: " << name << endl;
    break;
  case vault_severe:
  default:
    cout << " +++ DataVault Error:   Did not find item: " << name << endl;
    error("Fatal Error. ");
  }
  return -1;
}

String  DataVault::retrieve (const String& name,DataVaultCode code)
{
  int pos;
  if ((pos = find(name,v_string,code)) >= 0)
    return ncont[pos];
  return String("");
}

int  DataVault::insert(const String& name,const String& item,
		       DataVaultType type,DataVaultCode code)
{
  int pos = find(name);
  
  if (pos >= 0)
  {
    switch(code)
    {
    case vault_ignore:
      break;
    case vault_warning:
      cout << " +++ DataVault Warning: Item: " << name 
	   << " already present. " << endl;
      break;
    case vault_severe:
    default:
      cout << " +++ DataVault Error: Item: " << name 
	   << " already present." << endl;
      error("Fatal Error. ");
    }
    if (ntype[pos] != type)
    {
      cout << " +++ DataVault Error: Illegal to overwrite type for " << name 
	   << " Type was: " << ntype[pos] << " New Type: " << type << endl;
      error("Fatal Error. ");
    }
    ncont[pos] = item;
  }
  else
  {
    pos = nname.add(name);
    ntype.add(type);
    ncont.add(item);
  }
  return pos;
}

void   DataVault::unpack(char** ptr,double& d)
{
  char* str = *ptr;
  d = strtod(str,ptr);
  if (*ptr == str)
  {
    cout << "DataVault::Error Unpacking to double : " << str << endl;
    error("DataVault::Read Error. ");
  }
}

void      DataVault::unpack(char** ptr,long& d)
{
  char* str = *ptr;
  d = strtol(str,ptr,0);
  if (*ptr == str)
  {
    cout << "DataVault::Error Unpacking to long : " << str << endl;
    error("DataVault::Read Error. ");
  }
}

void    DataVault::pack(String& string,double d)
{
  char buf[100];
  sprintf(buf," %20.15e",d);
  String news(buf);
  string = string + news;
}

void    DataVault::pack(String& string,long d)
{
  char buf[100];
  sprintf(buf," %20li ",d);
  String news(buf);
  string = string + news;
}

void    DataVault::insert(const String& name,double  d,DataVaultCode code)
{
  String s;
  pack(s,d);
  insert(name,s,v_double,code);
}

void    DataVault::retrieve  (const String& name,double& d)
{
  int pos = find(name,v_double,vault_severe);
  char* item = ncont[pos].contents();
  char *start = item;
  unpack(&item,d);
  delete start;
}

void    DataVault::insert(const String& name,long d,DataVaultCode code)
{
  String s;
  pack(s,d);
  insert(name,s,v_int,code);
}

void    DataVault::retrieve(const String& name,int& d)
{
  int pos = find(name,v_int,vault_severe);
  char* item = ncont[pos].contents();
  char *start = item;
  long l;
  unpack(&item,l);
  d = (int) l;
  delete start;
}

void    DataVault::retrieve(const String& name,long& d)
{
  int pos = find(name,v_int,vault_severe);
  char *item  = ncont[pos].contents();
  char *start = item;
  unpack(&item,d);
  delete start;
}

void    DataVault::insert(const String& name,const NUMARRAY<int>& vec,
			  DataVaultCode code)
{
  String s;
  pack(s,(long) vec.size());
  
  for(int i =0; i < vec.size(); i++)
    pack(s,(long) vec(i));

  insert(name,s,v_ivector,code);
}

void    DataVault::retrieve(const String& name,NUMARRAY<int>& vec)
{
  int pos = find(name,v_ivector,vault_severe);
  char* item  = ncont[pos].contents();
  char* start = item;
  long no_item,v;
  unpack(&item,no_item);
  vec.reset(no_item);
  for(int i = 0; i < no_item; i++)
  {  
    unpack(&item,v);
    vec[i] = (int) v;
  }
  delete start;
}

void    DataVault::insert(const String& name,const Vec& vec,
			  DataVaultCode code)
{
  String s;
  pack(s,(long) vec.size());
  for(int i =0; i < vec.size(); i++)
    pack(s,vec(i));
  insert(name,s,v_vector,code);
}

void    DataVault::retrieve(const String& name,Vec& vec)
{
  int pos = find(name,v_vector,vault_severe);
  char* item  = ncont[pos].contents();
  char* start = item;
  long no_item;
  double v;
  unpack(&item,no_item);
  vec.reset(no_item);
  for(int i = 0; i < no_item; i++)
  {  
    unpack(&item,v);
    vec[i] = v;
  }
  delete start;
}


void    DataVault::insert(const String& name,const NUMARRAY2D<int>& mat,
			  DataVaultCode code)
{
  String s;
  pack(s,(long) mat.size1());
  pack(s,(long) mat.size2());
  for(int i =0; i < mat.size1(); i++)
  for(int j =0; j < mat.size2(); j++)
    pack(s,(long) mat.get(i,j));
  insert(name,s,v_imatrix,code);
}

void    DataVault::retrieve(const String& name,NUMARRAY2D<int>& mat)
{
  int pos = find(name,v_imatrix,vault_severe);
  char* item  = ncont[pos].contents();
  char* start = item;
  long d1,d2,v;
  unpack(&item,d1);
  unpack(&item,d2);
  mat.reset(d1,d2);
  for(int i = 0; i < d1; i++)
  for(int j = 0; j < d2; j++)
  {  
    unpack(&item,v);
    mat(i,j) = (int) v;
  }
  delete start;
}


void    DataVault::insert(const String& name,const Mat& mat,
			  DataVaultCode code)
{
  String s;
  pack(s,(long) mat.size1());
  pack(s,(long) mat.size2());
  for(int i =0; i < mat.size1(); i++)
  for(int j =0; j < mat.size2(); j++)
    pack(s,(double) mat.get(i,j));
  insert(name,s,v_matrix,code);
}

void    DataVault::retrieve(const String& name,Mat& mat)
{
  int pos = find(name,v_matrix,vault_severe);
  char* item  = ncont[pos].contents();
  char* start = item;
  long d1,d2;
  double v;
  unpack(&item,d1);
  unpack(&item,d2);
  mat.reset(d1,d2);
  for(int i = 0; i < d1; i++)
  for(int j = 0; j < d2; j++)
  {  
    unpack(&item,v);
    mat(i,j) = v;
  }
  delete start;
}



ostream& operator<<( ostream& os, DataVaultType& typ )
{

  switch(typ)
  {
  case v_string : os << "s "; break;
  case v_boolean: os << "b "; break;
  case v_int    : os << "i "; break;
  case v_double : os << "d "; break;
  case v_vector : os << "v "; break;
  case v_ivector: os << "p "; break;
  case v_matrix : os << "m "; break;
  case v_imatrix: os << "q "; break;
  default:
    error("DataVault::illegal type.");
  }
  return os;
  
}
  

istream& operator>>(istream& is,DataVaultType& t)
{
  char c;
  is >> c;
  if (!is)
    error("DataVault::Could not read type.");
  
  switch(c)
  {
  case 's': t =   v_string; break;
  case 'b': t =  v_boolean; break;
  case 'i': t =      v_int; break;
  case 'd': t =   v_double; break;
  case 'v': t =   v_vector; break;
  case 'p': t =  v_ivector; break;
  case 'm': t =   v_matrix; break;
  case 'q': t =  v_imatrix; break;
  case 'a': t =  v_svector; break;
  default:
    cout << " +++ DataVault:: read " << c << endl;
    error("DataVault::illegal type.");
  }
  return is;
}


void    DataVault::insert(const String& name,const ARRAY<String>& vec,
			  DataVaultCode code)
{
  String s;
  pack(s,(long) vec.size());
  for(int i =0; i < vec.size(); i++)
  {
    pack(s,(long) vec(i).length());
    s = s + vec(i);
  }  
  insert(name,s,v_svector,code);
}

void    DataVault::retrieve(const String& name,ARRAY<String>& vec)
{
  int pos = find(name,v_svector,vault_severe);
  char* item  = ncont[pos].contents();
  char* start = item;
  long no_item,v,len;

  unpack(&item,no_item);
  for(int i = 0; i < no_item; i++)
  {  
    unpack(&item,len);
    char tmp  = item[len];
    item[len] = 0;
    vec[i] = String(item);
    item[len] = tmp;
    item += len;
  }
  delete start;
}
