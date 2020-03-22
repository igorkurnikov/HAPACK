/*! \file hamatdb.cpp

    Classes to define Database of operator submatricies
 
    \author Igor Kurnikov 
    \date 1998-2003

*/

#define HAMATDB_CPP

#include <boost/algorithm/string.hpp>

#include "hamatdb.h"
#include "tokens.h"

   
StrKey::StrKey()
{

}

StrKey::StrKey(const char* str)
{
	SetStr(str);
}

StrKey::~StrKey()
{

}

void StrKey::SetStr(const char* str)
{
  key_str = str;
  boost::trim(key_str);
}


HaGrpOperID::HaGrpOperID()
{

}

HaGrpOperID::HaGrpOperID(const char* oper_type)
{
	this->set(oper_type);
}



HaGrpOperID::HaGrpOperID(const char* oper_type, const ChemGroup& g1, const ChemGroup& g2)
{
	this->set(oper_type,g1,g2);
}

HaGrpOperID::~HaGrpOperID()
{

}

void HaGrpOperID::set(const char* oper_type)
{	
	std::string tmp = oper_type;
	tmp += " 999 999";
	SetStr(tmp.c_str());
}


void HaGrpOperID::set(const char* oper_type, const ChemGroup& g1, const ChemGroup& g2)
{
	std::string tmp = oper_type;
  tmp += " ";
  tmp += g1.GetID();
  tmp += " ";
  tmp += g2.GetID();
  SetStr(tmp.c_str()); 
}


ostream& operator<<(ostream& s, HaGrpOperID& g_op_id)
{
	return(g_op_id.Print_To(s));
}


HaGrp4MatID::HaGrp4MatID()
{

}

HaGrp4MatID::HaGrp4MatID( const HaGrpOperID& g_op_id1, const HaGrpOperID& g_op_id2)
{
	key_str = g_op_id1.GetStr();
	key_str += "_";
	key_str += g_op_id2.GetStr();
}

 
HaGrp4MatID::~HaGrp4MatID()
{

}


HaMatDB::HaMatDB()
{  
   open_flag = 0; 
   debug_level=5;
}

HaMatDB::HaMatDB(const char* fname,const char* mode)
{
	debug_level=5;
	this->open(fname,mode);
}


HaMatDB::~HaMatDB()
{
	if(open_flag) this->close();
}

bool HaMatDB::open(const char* fname,const char* mode )
{
    open_flag = 1;
	return true;
}

bool HaMatDB::close()
{
    open_flag = 0;
	return true;
}

bool HaMatDB::is_open()
{
  if(open_flag) return true;
  return false;
}


bool HaMatDB::put(StrKey& key, const HaMat_double& fgmat)
{
	ostrstream data_buf_stream;
	data_buf_stream << fgmat;

	return true;
}

bool HaMatDB::get(StrKey& key, HaMat_double& fgmat)
{
	std::string data_str;
	try
	{
		if( data_str.empty() )
		{
			if(debug_level > 10)
			{
				cerr << " HaMatDB::get(): Submatrix not found with key: " << endl;
				key.Print_To(cout) ;
			}
			return false;
		}
		else
		{
			istrstream data_buf_stream( data_str.c_str() );
			size_t M, N;
            data_buf_stream >> M >> N;
			
			if( M > 0 && N > 0)
			{
			   fgmat.newsize(M,N);
               size_t nn = M*N;
			   size_t i,j;
			   for(i=0; i< M; i++)
			   {
			      for( j = 0; j < N ; j++)
			      {
			         double val;
                     data_buf_stream >> val;
                     fgmat.r0(i,j) = val;
				  }
			   }
			}
//			data_buf_stream >> fgmat;
		}
	}
	catch( const std::exception& ex)
	{
		cerr << ex.what() << "\n";
		return false;
	}
	
	return true; 	
}

bool HaMatDB::ListKeys() 
{
	return true;
 
}

bool HaMatDB::ListAll(ostream& sout) 
{
	return true;
 
}

int HaMatDB::PutMat(const char* key_str, const HaMat_double& fmat)
{
	StrKey key(key_str);
	return (int) put(key,fmat);
}

int HaMatDB::GetMat(const char* key_str, HaMat_double& fmat)
{
	StrKey key(key_str);
	return (int) get(key,fmat);
}

int  HaMatDB::ExtractToFile(const StrVec& keys, const char* fname)
{
	FILE* fout = fopen(fname,"w");
	if(fout == NULL) return FALSE;
	int nk = keys.size();
	int ik;
	HaMat_double dmat;
	for(ik = 0; ik < nk; ik++)
	{
		std::string key_str = keys[ik];
		GetMat(key_str.c_str(),dmat);
		int nr = dmat.num_rows();
		int nc = dmat.num_cols();
		if( nr > 0 && nc > 0)
		{
			fprintf(fout,"%s\n",key_str.c_str());
			fprintf(fout,"%6d %6d \n", nr,nc);
			int i,j;
			for(i = 0; i < nr; i++)
			{
				for(j = 0; j < nc; j++)
				{
					fprintf(fout," %12.5f ",dmat.GetVal_idx0(i,j));
				}
				fprintf(fout,"\n");
			}
		}
	}
	fclose(fout);
	return TRUE;
}

int HaMatDB::ExtractAllToFile(const char* fname)
{
	if(!this->is_open()) return false;
	    

    // Walk through the Database, printing the keys.
    StrKey key;
	StrVec keys;
 
	return ExtractToFile(keys,fname);
}

int 
HaMatDB::AddFromFile(const char* fname)
{
	char buf[2400];

	ifstream is(fname);
	if( is.fail() )
	{	
		PrintLog(" Error in HaMatDblDB::AddFromFile() \n");
		PrintLog(" Error opening file %s \n",fname);
		return FALSE;
	}

	for(;;)
	{
		is.getline(buf,255);
		if( is.fail() ) return TRUE;
		std::string key_str = buf;
		HaMat_double fmat;
		is.getline(buf,255);
		if(is.fail()) 
		{
			PrintLog(" Error reading matrix for key %s \n",key_str.c_str());
			return FALSE;
		}
		istrstream is2(buf);
		int nr, nc;
		is2 >> nr;
		is2 >> nc;
		PrintLog(" nr = %d  nc = %d \n",nr,nc);
		if( is2.fail())
		{
			PrintLog(" Error reading matrix dimensions for key %s",key_str.c_str());
			return FALSE;
		}
		int i,j;
		fmat.newsize(nr,nc);
		for( i=0; i < nr; i++)
		{
			is.getline(buf,255);
			if( is.fail() )
			{
				PrintLog(" Error reading matrix values for key %s",key_str.c_str());
				return false;
			}
			istrstream is3(buf);
			for(j = 0; j < nc; j++)
			{
				double val;
				is3 >> val;
		    	if( is3.fail() )
				{
				   PrintLog(" Error reading matrix values for key %s",key_str.c_str());
				   return false;
				}
				fmat.SetVal_idx0(i,j,val);
			}
		}
		PutMat(key_str.c_str(),fmat);
	}
	return TRUE;
}



