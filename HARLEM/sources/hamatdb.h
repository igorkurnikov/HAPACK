/*! \file hamatdb.h

    Classes to define a database of operator submatricies in HARLEM

   \author Igor Kurnikov
   \date 1998-2008

*/
#ifndef HAMATDB_H
#define HAMATDB_H

#include "halinalg.h"
#include "haatgroup.h"
#include "command.h"

class StrKey
{
public:
	StrKey();
	StrKey(const char* str);
	virtual ~StrKey();

	void SetStr(const char* str);
	const char* GetStr() const { return key_str.c_str();} 

	std::ostream& Print_To( std::ostream& os ) { return (os << key_str); }

	std::string key_str;
};


class HaGrpOperID: public StrKey
{
public:
	HaGrpOperID();
	HaGrpOperID(const char* oper_type);
	HaGrpOperID(const char* oper_type, 
		        const ChemGroup& g1, const ChemGroup& g2);
	virtual ~HaGrpOperID();

    void set(const char* oper_type);
	void set(const char* oper_type, 
    const ChemGroup& g1, const ChemGroup& g2);

};

class HaGrp4MatID: public StrKey
//!< class for ID of 4-group index matrix 
{
public:
	HaGrp4MatID();
	HaGrp4MatID( const HaGrpOperID& g_op_id1, const HaGrpOperID& g_op_id2);
	virtual ~HaGrp4MatID();

};

class HaMatDB 
//!  Class to define Database of sub matricies
{
public:
  HaMatDB();
  HaMatDB(const char* fname,const char* mode="r" );
  virtual ~HaMatDB();

  virtual bool open(const char* fname, const char* mode );
  bool close();

  bool is_open(); 

  bool put(StrKey& key, const HaMat_double& fgmat);  
  bool get(StrKey& key, HaMat_double& fgmat);
  
  virtual bool ListKeys() ;
  virtual bool ListAll(ostream& s);

  int PutMat(const char* key_str, const HaMat_double& fmat );
  int GetMat(const char* key_str, HaMat_double& fmat);
	
  int ExtractToFile(const StrVec& keys, const char* fname);
  int ExtractAllToFile(const char* fname);
  int AddFromFile(const char* fname);

  std::string io_format;

protected:
  int open_flag;
  int debug_level;
  
  HaVec_int data_buffer;
};

#endif /* !HAMATDB_H */

