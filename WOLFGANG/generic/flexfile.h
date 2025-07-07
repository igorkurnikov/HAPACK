/*TEX
\section{Flexible File Structure}

Associated with a \name{FlexFile} are:

an index-file, which for each key gives a list of buckets (name: basename.idx)


*/

#if !defined(_MSC_VER)
#include <unistd.h>
#endif
#include <tlist.h>
#include <tsstack.h>
#include "corefile.h"
#include "io_incl.h"

template <class Value> class FlexFile : virtual public BucketFile<Value>
{
  char*                          basename;
  int                            current_item;
  size_t                         bsize;         // size of bucket field in char
  FILE*                          file;

  SList < int >                  freed_items;
  ARRAY < SList  <int> >         index;

  Bucket<Value>                  read_buffer;

void         write_index();
void         read_index();

public:
          /// open existing FlexFile 
          FlexFile(const char* basename,int debug_mode = 0);  
          /// open new FlexFile
          FlexFile(const char* basename,int no_keys,
		   int bucket_size,int debug_mode = 0);
          /// close FlexFile
virtual  ~FlexFile() { finalize(); };
	  /// add an item
        void write_all_basic();
        void write_bucket(const int k);
virtual void compress();  
virtual int  local(int) { return 1; }  
virtual int  addr (int) { return 0; }  
virtual const Bucket<Value>& read(int key,int unlink = 0);
void    finalize();
};


template <class Value> FlexFile<Value>::FlexFile(const char* basen,
						 int no_keys,int b_size,
						 int debug_mode)
  : BucketFile<Value>(no_keys,b_size,debug_mode)
{
  current_item  = 0;
  bsize         = 2*sizeof(int) + sizeof(Value)*this->bucket_size;
  basename      = strdup(basen);
  
  // this is a new file
  
  index .reset (no_keys);
  char buf[100];
  sprintf(buf,"%s.dat",basen);
  cout << " +++ Opening FlexFile named: " << buf << endl;
  
  file = fopen(buf,"w+");
  if (!file)
    error(" +++ FlexFile<Value>::Error =-= Could not open FlexFile ");
}

template <class Value> FlexFile<Value>::FlexFile(const char* base,
						 int debug_mode)
  : BucketFile<Value>(0,0,debug_mode)
{
  basename = strdup(base);  
  char buf[100];
  sprintf(buf,"%s.dat",basename);
  cout << " +++ Opening FlexFile named: " << buf << endl;
  
  file = fopen(buf,"r+");
  if (!file)
    error(" +++ FlexFile<Value>::Error =-= Could not open FlexFile ");

  read_index();
}

/*TEX
\subsection{Writing the index file}

Writing the index file is a destructive operation. The old contents is wholesale replaced.
Writing the index file entails a flushing of the buffers. 

*/

template <class Value> void FlexFile<Value>::write_index()
{
  char buf[200];
  sprintf(buf,"%s.idx",basename);
  FILE* index_file = fopen(buf,"w");
  if (!index_file)
    error("FlexFile<Value>::Could not open index file.");

  int no_keys = index.size();
  int code = fwrite(&no_keys,sizeof(int),1,index_file);
  if (code != 1)
    error("FlexFile<Value>::Could not write index file.");

  code = fwrite(&this->bucket_size,sizeof(int),1,index_file);
  if (code != 1)
    error("FlexFile<Value>::Could not write index file.");

  code = fwrite(&current_item,sizeof(int),1,index_file);
  if (code != 1)
    error("FlexFile<Value>::Could not write index file.");

  for(int i = 0; i < index.size(); i++)
    index[i].write(index_file);

  fclose(index_file);
}

template <class Value> void FlexFile<Value>::read_index()
{
  char buf[200];
  sprintf(buf,"%s.idx",basename);
  FILE* index_file = fopen(buf,"r");
  if (!index_file)
    error("FlexFile<Value>::Could not open index file.");

  int no_keys;
  int code = fread(&no_keys,sizeof(int),1,index_file);
  if (code != 1)
    error("FlexFile<Value>::Could not write index file.");

  code = fread(&this->bucket_size,sizeof(int),1,index_file);
  if (code != 1)
    error("FlexFile<Value>::Could not write index file.");

  code = fread(&current_item,sizeof(int),1,index_file);
  if (code != 1)
    error("FlexFile<Value>::Could not write index file.");

  bsize         = 2*sizeof(int) + sizeof(Value)*this->bucket_size;

  index .reset(no_keys);
  this->bucket.reset(no_keys);

  for(int i = 0; i < index.size(); i++)
    index[i].read(index_file);

  fclose(index_file);
}

/*TEX 
\subsection{Flushing the Buffers}
*/
template <class Value> void FlexFile<Value>::write_all_basic()
{
  FlexFile<Value>* ff = this;
  for(int i = 0; i < this->bucket.size(); i++)
  {
    int code;
    if (this->bucket[i].size() > 0)
    {
      long pos; 
      
      if (freed_items.size() > 0) // use freed item
      {
	pos = freed_items();
	freed_items.remove();
      }
      else 
	pos = current_item++;
      
      index[i].insert(pos);
      
      if (this->_debug)
	cout << " Bucket for Key: " << setw(3) << i 
	     << " pos: " << setw(5) << pos 
	     << " sz: " << index[i].size() << endl;
      
      code = fseek(file,bsize * pos,SEEK_SET);
      if (code != 0)
      {
	cout << " +++ FlexFile::Error =-= seeking bucket " 
	     << pos << endl;
	error(" +++ FlexFile::Error seeking bucket in file.");
      }
      
      code = fwrite(&i,sizeof(int),1,file);
      if (code != 1)
      {
	cout << " +++ FlexFile::Error =-= writing bucket id " 
	     << pos << endl;
	error(" +++ FlexFile::Error writing bucket index to file. ");
      }
      
      int cno = this->bucket[i].size();
      code = fwrite(&cno,sizeof(int),1,file);
      if (code != 1)
      {
	cout << " +++ FlexFile::Error =-= writing bucket id " 
	     << pos << endl;
	error(" +++ FlexFile::Error writing bucket index to file. ");
      }
      
      code = fwrite(this->bucket[i].get_base(),sizeof(Value),this->bucket_size,file);
      if (code != this->bucket_size)
      {
	cout << " +++ FlexFile::Error =-= writing bucket " 
	     << pos << endl;
	error(" +++ FlexFile::Error writing bucket to file. ");
      }
      this->bucket[i].reset();
    }
  }
}


/*TEX
\subsection{Writing a bucket}

We determine a free position, either from the freed\_items list or at the end of the file.
We then write the key, the actual size and the entire bucket, independent of its actual filling.
*/

template <class Value> void FlexFile<Value>::write_bucket(const int i)
{
  int code;
  if (this->bucket[i].size() > 0)
  {
    long pos; 

    if (freed_items.size() > 0) // use freed item
    {
      pos = freed_items();
      freed_items.remove();
    }
    else 
      pos = current_item++;

    index[i].insert(pos);

    if (this->_debug)
      cout << " Bucket for Key: " << setw(3) << i << " pos: " << setw(5) << pos 
	   << " sz: " << index[i].size() << endl;
 
    code = fseek(file,bsize * pos,SEEK_SET);
    if (code != 0)
    {
      cout << " +++ FlexFile::Error =-= seeking bucket " << pos << endl;
      error(" +++ FlexFile::Error seeking bucket in file.");
    }

    code = fwrite(&i,sizeof(int),1,file);
    if (code != 1)
    {
      cout << " +++ FlexFile::Error =-= writing bucket id " << pos << endl;
      error(" +++ FlexFile::Error writing bucket index to file. ");
    }
 
    int cno = this->bucket[i].size();
    code = fwrite(&cno,sizeof(int),1,file);
    if (code != 1)
    {
      cout << " +++ FlexFile::Error =-= writing bucket id " << pos << endl;
      error(" +++ FlexFile::Error writing bucket index to file. ");
    }
 
    code = fwrite(this->bucket[i].get_base(),sizeof(Value),this->bucket_size,file);
    if (code != this->bucket_size)
    {
      cout << " +++ FlexFile::Error =-= writing bucket " << pos << endl;
      error(" +++ FlexFile::Error writing bucket to file. ");
    }
    this->bucket[i].reset();
  }
}

/*TEX
\subsection{Compress}

Compression is a non-trivial operation, as the index does not yield the information in the 
desired format. We must move buckets from the end of the file to freed positions, such that

 (i) all freed positions, which are followed by at least one nonempty bucket are filled and
(ii) only freed positions at the end of the file remain. 

This is accomplished in the following steps: 

(i)   count the number of exisiting items. All positions from 0 to no\_item-1 
      have then to be filled. All freed positions with larger indices can be discared.
(ii)  run through the list of all buckets. If a bucket is in a position larger or equal 
      than no\_item, copy it to the next free positions less than no\_item. 
(iii) repeat the previous steps until there are no more free positions.
(iv)  truncate the file. 

*/

template <class Value> void FlexFile<Value>::compress()
{
  write_all_basic(); // requires all buckets to be flushed.
  int no_items = 0;
  int k;
  for(k = 0; k < index.size(); k++)
    no_items += index[k].size();

  char* buf = new char[bsize];
  if (!buf)
    error("ENOMEM in compress.");
  
  if (this->_debug)
    cout << " +++ Compress: Number of Items: " << no_items << endl;
  
  int code;

  while (freed_items.size() > 0 && freed_items() >= no_items)
    freed_items.remove();

  for(k = 0; k < index.size(); k++)
  {
    if (freed_items.size() == 0) break; // nothing left to do
    for(SListIterator<int> iter(index[k]);iter.in_range(); ++iter)
    {
      if (freed_items.size() == 0) break; // nothing left to do
      if (iter() >= no_items) // element to write
      {   
	int pos = freed_items();

	if (this->_debug)
	  cout << " --- compress moving bucket with key: " << k << " from position: " 
	       << iter() << " to position: " << pos << endl;
 
	freed_items.remove();
	while (freed_items.size() > 0 && freed_items() >= no_items)
	  freed_items.remove();
  
	code = fseek(file,bsize *  iter(),SEEK_SET);
	if (code != 0)
	{
	  cout << " +++ FlexReader::Error =-= seeking bucket R: " << iter() 
	       << " in compress. " << endl;
	  error(" +++ FlexFile::Error seeking bucket in file.");
	}
	code = fread(buf,sizeof(char),bsize,file);
	if (code != bsize)
	{
    	  cout << " +++ FlexReader::Error =-= reading bucket " << iter() 
	       << " in compress. " << endl;
	  error(" +++ FlexFile::Error seeking bucket in file");
	}

	code = fseek(file,bsize *  pos,SEEK_SET);
	if (code != 0)
	{
	  cout << " +++ FlexReader::Error =-= seeking bucket W: " << pos
	       << " in compress. " << endl;
	  error(" +++ FlexFile::Error seeking bucket in file");
	}

	code = fwrite(buf,sizeof(char),bsize,file);
	if (code != bsize)
	{
    	  cout << " +++ FlexReader::Error =-= writing bucket " << pos
	       << " in compress. " << endl;
	  error(" +++ FlexFile::Error seeking bucket in file");
	}
	iter() = pos; // remember the new position
      }
    }
  }
  
  fclose(file);
  char fname[100];
  sprintf(fname,"%s.dat",basename);
 
#if !defined(_MSC_VER)
  code = truncate(fname,no_items*bsize);
  if (code != 0)
    error("FlexFile =-= Truncating File. ");
#endif
  current_item = no_items;

  file = fopen(fname,"r+");
  if (!file)
    error("FlexFile =-= could not reopen file in compress.");

  delete buf;
}

/*TEX
\subsection{Destructor}
*/
template <class Value> void FlexFile<Value>::finalize()
{
  cout << " +++ Closing Integral File. " << endl;
  write_all_basic(); 
  write_index();
  fclose(file);
  freed_items.clear();
  int i;
  for(i = 0; i < index.size(); i++)
    index[i].clear();
  for(i = 0; i < this->bucket.size(); i++)
    this->bucket[i].resize(0);
  
  delete basename;
}

template <class Value> const Bucket<Value>& FlexFile<Value>::read(int key,
								  int unlink) 
// third argument is ignored
{
  // count the number of items
  
  int mx_sz = this->bucket_size * index[key].size() + this->bucket[key].size();
  if (read_buffer.max_size() < mx_sz)
    read_buffer.resize(mx_sz);
  read_buffer.reset();
  
  int sz = 0; // keep track of items actually read.
  int code;

  for(SListIterator<int> iter(index[key]); iter.in_range(); ++iter)
  {

    code = fseek(file,bsize * iter(),SEEK_SET);
    if (code != 0)
    {
      cout << " +++ FlexReader::Error =-= seeking bucket for key: " << key 
	   << " in position " << iter() << endl;
      error(" +++ FlexFile::Error seeking bucket in file.");
    }

    int kk;
    code = fread(&kk,sizeof(int),1,file);
    if (code != 1)
    {
      cout << " +++ FlexFile::Error =-= reading bucket for key: " << key  
	   << " in position " << iter() << endl;
      ::error(" +++ FlexFile::Error reading bucket index from file.");
    }
        
    if (kk != key)
    {
      cout << " +++ FlexFile::Error =-= bucket key mismatch for bucket no: " << iter() 
	   << " had: " << kk << " expected: " << key << endl;
      // error(" +++ FlexFile::Error bucket key mismatch.");
    }
      
    int bsize;
    code = fread(&bsize,sizeof(int),1,file);
    if (code != 1)
    {
      cout << " +++ FlexFile::Error =-= reading bucket size for key " << key 
	   << " in position: " << iter() << endl;
      error(" +++ FlexFile::Error reading bucket index to file.");
    }
      
    code = fread(read_buffer.get_base() + sz,sizeof(Value),bsize,file);
    if (code != bsize)
    {
      cout << " +++ FlexFile::Error =-= reading bucket for key: " << key 
	   << " in position " << iter() << endl;
      ::error(" +++ FlexFile::Error reading bucket.");
    }
    sz += bsize;

    if (unlink)
      freed_items.insert(iter()); 
  }
    
  SStack<Value,0>& stack  = this->bucket[key];
// No function set_item_number() !! Incompatible lib??
//  read_buffer.set_item_number(sz);
  read_buffer.resize(sz);
  for(int i = 0; i < stack.size(); i++)
    read_buffer.add(stack[i]);
  
  if (unlink)
  {
    stack.resize(0);
    index[key].clear();
  }

  return read_buffer;
  
}
