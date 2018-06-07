/*TEX
\section{Distributed Core File}

if the ParMachine pointer pm is zero, a single-node version is assumed
*/
#ifndef COREFILE_H
#define COREFILE_H
#if !defined(_MSC_VER)
#include <unistd.h>
#endif
#include <tlist.h>
#include <tsstack.h>
#include "parallel.h"
#include "io_incl.h"


template <class Value> class Bucket : public SStack<Value,0>
{
public:
    Bucket()  {}
   ~Bucket()  {}
};
/*TEX

*/
template<class Value> ostream& operator<<(ostream& os, 
					  const Bucket<Value>& b)
{
  os << "too dumb to print Bucket";
  return os;
}

template <class Value> class BucketFile 
{
protected:
  int                            bucket_size;
  ARRAY <Bucket<Value> >         bucket; 
  int                            _debug;
  
public:
          BucketFile(int no_keys,int bucket_size,
		     int debug_mode);
virtual  ~BucketFile();
void      add(const int key,const Value& val);

virtual const Bucket<Value>& read (int key,int unlink) = 0;
virtual void          write_bucket(const int key)      = 0;
        void          clean_buckets(); 
};
 
template <class Value> BucketFile<Value>::BucketFile(int no_keys,int b_size,
						     int debug_mode)
{
  bucket_size = b_size;  
  bucket.reset(no_keys);
  _debug = debug_mode;
}

template <class Value> void BucketFile<Value>::add(const int key,
						   const Value& val)
{
  Bucket<Value>& bb = bucket[key];
  if (bb.max_size() == 0)
    bb.resize(bucket_size);
  bb.add(val);
  if (bb.size() == bucket_size)
    write_bucket(key);
}

template <class Value> BucketFile<Value>::~BucketFile()
{
  clean_buckets();
}

template <class Value> void BucketFile<Value>::clean_buckets()
{
  int i;
  for(i = 0; i < bucket.size(); i++)
  {
    if (bucket[i].size() > 0)
      cout << " +++ Warning: Bucket "  << i << " is not empty " << endl;
    bucket[i].resize(0);  
  }
  
}
/*TEX
This class must implement the virtual functions in BucketFileHandle
*/
template <class Value> class CoreFile : virtual public BucketFile<Value>
{
  ARRAY <SList<ARRAY<Value>*> >  local_store;  // local store
  ARRAY <SList<ARRAY<Value>*> >  global; // global store
  int                            rcvd;
  IVec                           address;
  Bucket<Value>                  read_buffer;
  ARRAY <Value>                  buffer;
  
public:
  CoreFile(int no_keys,int b_size,int debug_mode = 0);
  CoreFile(char* fname);
  ~CoreFile();
  void                         write_bucket(const int key);
  virtual const Bucket<Value>& read(int key,int unlink= 0);
  void                         write_all_basic(); 
  virtual void                 compress()     { write_all_basic();}
  virtual int                  local(int key) { return address(key) == pm -> id(); }
  virtual int                  addr (int key) { return address(key); }
  virtual int                  number_of_keys(){return address.size();}
  
  void                         write_file(const char* fname,int clear);
  void                         read_file (const char* fname);  
  void                         clean();  
};

template <class Value> CoreFile<Value>::CoreFile(int no_keys,int b_size,
						 int debug_mode)
  : BucketFile<Value>(no_keys,b_size,debug_mode)
{
  if (!pm)
    error("Global Pointer ParMachine::pm must be set to instantiate CoreFile");
  
  this->_debug        = debug_mode;
  buffer.reset(this->bucket_size);
  local_store.reset(no_keys);
  global.reset(no_keys);
  address.reset(no_keys);

  if ((this->_debug & 2) && pm -> master())
    cout << "  +++ Address Table: " << endl;

  /*
    int ipn  = no_keys/pm->nodes() + 1; // integs per node
    from_key = ipn*pm->id();
    to_key   = ipn*(pm->id()+1);
    if (to > no_keys) to = no_keys;
  */
  for(int i = 0; i < no_keys; i++)
  {
    address[i] = i % pm -> nodes();
    if ((this->_debug & 2) && pm -> master())
      cout << "        " << setw(4) << i << setw(4) << address(i) << endl;
  }  
  rcvd = 0;
}
/*TEX 
\subsection{Local File IO}

write all local data
*/
template <class Value> void CoreFile<Value>::write_file(const char* fname,
						   int clear)
{
  int code;
  char buf[100];
  FILE* file = fopen(fname,"w");
  if (!file)
    error("Opening Outputfile for Corefile");
  
  for(int key = 0; key < global.size(); key++)
  {
    const Bucket<Value>& reader = read(key);
    if (reader.size() > 0)
    {
      code = fwrite(&key,sizeof(int),1,file);
      if (code != 1)
	error("write error for key in corefile !");
      int sz = reader.size();
      code = fwrite(&sz,sizeof(int),1,file);
      if (code != 1)
	error("write error for sz in corefile !");
      code = fwrite(reader.get_base(),sizeof(Value),sz,file);
      if (code != sz)
	error("write error for sz in corefile !");
    }
  }
}

template <class Value> void CoreFile<Value>::read_file(const char* fname)
{
  char buf[100];
  sprintf(buf,"%s%03i",fname,pm -> id());
  FILE* file = fopen(fname,"r");
  if (!file)
    error("Could not open file in Corefile read");
  int code;
  
  for(int key = 0; key < global.size(); key++)
  {
    const Bucket<Value>& reader = read(key);
    if (reader.size() > 0)
    {
      code = fwrite(&key,sizeof(int),1,file);
      if (code != 1)
	error("write error for key in corefile !");
      int sz = reader.size();
      code = fwrite(&sz,sizeof(int),1,file);
      if (code != 1)
	error("write error for sz in corefile !");
      code = fwrite(reader.get_base(),sizeof(Value),sz,file);
      if (code != sz)
	error("write error for sz in corefile !");
    }
  }
}

/*TEX
Cleaing the COREFILE of all Data
*/
template <class Value> void CoreFile<Value>::clean()
{
  cout << " +++ Cleaning IFile: " << endl;

  this->clean_buckets();
  int i;
  for(i = 0; i < local_store.size(); i++)
    if (local_store[i].size() > 0)
      cout << " Warning: Local Store elements left for key: " 
	   << i << endl;

  for(i = 0; i < global.size(); i++)
  {
    read(i,1); // remove
  }
  
  
}


/*TEX 
\subsection{Flushing the Buffers}

write all remaining data
*/

template <class Value> void CoreFile<Value>::write_all_basic()
{
  int i,id;

  IVec sendcount (pm -> nodes());
  IVec recvcount (pm -> nodes());
  IVec senddispl (pm -> nodes());
  IVec recvdispl (pm -> nodes());
  IVec recvdispl1(pm -> nodes());
  IVec all_items(local_store.size());

  // step 1 communicate to each node the number of items for
  // each bucket from each of the other nodes

  sendcount.set(0);
  for(id = 0; id < local_store.size(); id++)
  {
    write_bucket(id);
    sendcount[address[id]]++;
  }
  this->clean_buckets();
  
  int send_total = 0;
  for(i = 0; i < pm -> nodes(); i++)
  {
    senddispl[i] = send_total;
    send_total  += sendcount[i];
  }

  pm -> all_to_all(sendcount,recvcount,1);
  
  // compute the receive counts and offsets 

  int recv_total = 0;
  for(i = 0; i < pm -> nodes(); i++)
  {
    recvdispl1[i] = recv_total;
    recv_total  += recvcount[i];
  }

  // fill the send buffers

  IVec sendbuf1(send_total);
  IVec sendpos(senddispl);
  for(id = 0; id < local_store.size(); id++)
  {
    int sz = 0;
    for(SListIterator<ARRAY<Value>*> counter(local_store[id]); 
	counter.in_range(); ++counter)
      sz += counter() -> size();  
    if (this->_debug)
      cout << " From : " << setw(4) << pm -> id() << " for id " 
	   << setw(4) << id << " #items " << setw(5) << sz << endl;
    sendbuf1[sendpos[address[id]]++] = sz;
  }

  IVec recvbuf1(recv_total);
  pm -> all_to_all_displ(sendbuf1,sendcount,senddispl, 
			 recvbuf1,recvcount,recvdispl1);

  // unpack the send buffers
  int pos   = 0;
  int count = 0;
  all_items.set(0);
  for(id = 0; id < global.size(); id++)
    if (local(id))
    {
      if (this->_debug)
	cout << " Node: " << pm -> id() << " receiving for id: " 
	     << id << " : ";
      for(i = 0; i < pm -> nodes(); i++)
      {
	all_items[id] += recvbuf1[recvdispl1[i] + count];
	if (this->_debug)
	  cout << setw(4) << recvbuf1[recvdispl1[i] + count];
      }
      if (this->_debug)
	cout << " Tot: " << all_items[id] << endl;
      count++;
    }
  
  // step 2: transfer the data
  // compute the send counts and  offsets

  sendcount.set(0);
  for(id = 0; id < local_store.size(); id++)
  {
    int sz = 0;
    for(SListIterator<ARRAY<Value>*> counter(local_store[id]); 
	counter.in_range(); ++counter)
      sz += counter() -> size();  
    sendcount[address[id]] += sz;
  }
  
  send_total = 0;
  for(i = 0; i < pm -> nodes(); i++)
  {
    senddispl[i] = send_total;
    send_total  += sendcount[i];
    if (this->_debug)
      cout << " From : " << setw(4) << pm -> id() << " to: " 
	   << setw(4) << i << " #items " << setw(5) << sendcount[i]
	   << " displ: " << setw(5) << senddispl[i] << endl;
  }

  pm -> all_to_all(sendcount,recvcount,1);
  
  // compute the receive counts and offsets 

  recv_total = 0;
  for(i = 0; i < pm -> nodes(); i++)
  {
    recvdispl[i] = recv_total;
    recv_total  += recvcount[i];
    if (this->_debug)
      cout << " At : " << setw(4) << pm -> id() << " from: " 
	   << setw(4) << i << " #items " << setw(5) << recvcount[i] 
	   << " displ " << setw(5) << recvdispl[i] << endl;
  }
  if (this->_debug)
    for (id = 0; id < all_items.size(); id++)
      if (local(id))
	cout << " At : " << setw(4) << pm -> id() << " for: " 
	     << setw(4) << id << " #items " << setw(5) 
	     << all_items[id] << endl;

  // fill the send buffers

  ARRAY<Value> sendbuf(send_total);
  sendpos = senddispl;
  for(id = 0; id < local_store.size(); id++)
  {
    for(SListIterator<ARRAY<Value>*> iter(local_store[id]); 
	iter.in_range(); ++iter)
    {
      ARRAY<Value>& a = *(iter());
      for(int i = 0; i < a.size(); i++)
	sendbuf[sendpos[address[id]]++] = a[i];
      delete iter();
    }
    local_store[id].clear();
  }

  ARRAY<Value> recvbuf(recv_total);
  for(i = 0; i < pm -> nodes(); i++)
  {
    sendcount[i] *= sizeof(Value);
    recvcount[i] *= sizeof(Value);
    senddispl[i] *= sizeof(Value);
    recvdispl[i] *= sizeof(Value);
  }
  
  pm -> all_to_all_displ(sendbuf.get_base(),sendcount,senddispl,
			 recvbuf.get_base(),recvcount,recvdispl);
  
  recvdispl /= sizeof(Value);
  
  // unpack the send buffers

  pos   = 0;
  count = 0;
  sendpos = recvdispl;

  for(id = 0; id < global.size(); id++)
    if (local(id))
    {
      if (this->_debug)
	cout << " Node: " << pm -> id() << " receiving " << setw(5) 
	     << all_items[id] << " for id: " << setw(4) << id << endl;
      ARRAY<Value>* buf = new ARRAY<Value>(all_items[id]);
      pos = 0;      
      for(i = 0; i < pm -> nodes(); i++)
      {
	int items = recvbuf1[recvdispl1[i] + count];
	if (this->_debug)
	  cout << " .... " << pm -> id() << " / " << id << " from: " << i 
	       << " #items: " << items << " reading from: "
	       << sendpos[i] << " writing to: " << pos << endl;
	for(int k = 0; k < items; k++)
	  (*buf)[pos++] = recvbuf[sendpos[i]++];	  
      }
      if (pos != all_items[id])
	cout << " Node: " << pm -> id() << " expected: " << all_items[id] 
	     << " items for id: " << id << " but got: " << pos 
	     << " items. " << endl;
      count++; // counts the legal id's 
      global[id].insert(buf);
    }

  
  if (this->_debug)
    cout << " Node: " << pm -> id() << " finalized corefile. " << endl;
}

/*TEX
\subsection{Writing a bucket}

We determine a free position, either from the freed\_items list or at
the end of the file.  We then write the key, the actual size and the
entire bucket, independent of its actual filling.  
*/

template <class Value> void CoreFile<Value>::write_bucket(const int key)
{
  int code;
  if (this->bucket[key].size() > 0)
  {
    // local copy
    int no = this->bucket[key].size();
    ARRAY<Value>* data = new ARRAY<Value>(no);
    if (!data)
      error("ENOMEM for data in write-bucket of CoreFile");
    for(int i = 0; i < no; i++)
      (*data)[i] = this->bucket[key][i];      
    if (address(key) == pm -> id())
      global[key].insert(data);
    else 
      local_store[key].insert(data);
    this->bucket[key].reset();
  }
}

/*TEX
\subsection{Synchronize}
Receive a buffer and add data.
*/

/*template <class Value> void CoreFile<Value>::receive()
{
  int code,key;
  int no = bucket_size*sizeof(Value);
  code = pm -> receive((void*) buffer.get_base(),no);
  no  /= sizeof(Value);
  key  = pm -> tag();
  if (_debug)
    cout << " Node: " << pm -> id()
	 << " received: " << no << " for key: " << key
	 << " from: " << pm -> source() << endl;
  ARRAY<Value>* data = new ARRAY<Value>(no);
  if (!data)
    error("ENOMEM for data in CoreFile");  
  for(int i = 0; i < no; i++)
    (*data)[i] = buffer[i];
  local_store[key].insert(data);
  rcvd++;
}
*/  
/*TEX
\subsection{Destructor}
*/
template <class Value> CoreFile<Value>::~CoreFile()
{
  int i;
  for(i = 0; i < local_store.size(); i++)
  {
    while(local_store[i].size() > 0)
    {
      delete local_store[i]();
      local_store[i].remove();
    }
  }

  for(i = 0; i < global.size(); i++)
  {
    while(global[i].size() > 0)
    {
      delete global[i]();
      global[i].remove();
    }
  }
}


template <class Value> const Bucket<Value>& CoreFile<Value>::read(int key,int unlink)
{
  read_buffer.reset();
  
  if (address(key) != pm -> id()) // local data only
    return read_buffer;

  int sz = 0;
  for(SListIterator<ARRAY<Value>*> counter(global[key]); 
      counter.in_range(); ++counter)
    sz += counter() -> size();  
  if (read_buffer.max_size() < sz)
    read_buffer.resize(sz);

  int count = 0;
  for(SListIterator<ARRAY<Value>*> iter(global[key]); 
      iter.in_range(); ++iter)
  {
    ARRAY<Value>& a = *(iter());
    for(int i = 0; i < a.size(); i++)
      read_buffer.add(a[i]);
    if (unlink)
      delete iter();   
  }
  if (unlink)
    global[key].clear();

  return read_buffer;
}


#endif
