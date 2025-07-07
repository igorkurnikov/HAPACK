#include "parallel.h"
#include "io_incl.h"
#include <errorip.h>
#include <assert.h>
#include <tarray.h>
#include "vtype.h"
#include <numer.h>

ParMachine::ParMachine(int* arfc,char*** a) {}
ParMachine::~ParMachine() {}
int  ParMachine::all_sum  (IVec& src,IVec& tar) { tar = src; return 0; }
int  ParMachine::all_sum  (Mat & src,Mat & tar) { tar = src; return 0; }
int  ParMachine::all_sum  (NUMARRAY2D<float> & src,NUMARRAY2D<float> & tar) { tar = src; return 0; }
int  ParMachine::all_to_all(IVec& in,IVec& out,int items) {  out = in; return 0; }
int  ParMachine::all_to_all_displ(IVec&  in,IVec& sendcount,IVec& senddispl,
				  IVec& out,IVec& recvcount,IVec& recvdispl)
    { out = in; recvcount = sendcount; recvdispl = senddispl; return 0; }


ParMachine *pm;

#ifdef MPI
MPI_Group CompGroup;
MPI_Comm  CompDomain;

MPIParMachine::MPIParMachine(int* argc,char*** argv) : ParMachine(argc,argv)
{
  MPI_Init(argc,argv);
  my_id           = 0;
  numprocs        = 1;
  world_id        = 0;
  world_numprocs  = 1;
  
  MPI_Comm_size (MPI_COMM_WORLD,&world_numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&world_id);
  MPI_Comm_group(MPI_COMM_WORLD,&CompGroup);
#ifdef CSIM
  int rank = 0;
  MPI_Group_excl(CompGroup, 1, &rank,&CompGroup);
#else
  const int CSIM = 0;
#endif
  MPI_Comm_create(MPI_COMM_WORLD,CompGroup,&CompDomain);

  if (world_id || !CSIM)
  {
    MPI_Comm_rank(CompDomain,&my_id);
    MPI_Comm_size(CompDomain,&numprocs);
  }
  else 
  {
    numprocs = 0;
    my_id    = 0;
  }
  rq.reset(numprocs);
  rq.set(0);
}
MPIParMachine::~MPIParMachine()
{
  MPI_Finalize();
}

void MPIParMachine::abort(int code)
{
  MPI_Abort(MPI_COMM_WORLD,code);
}

int MPIParMachine::all_to_all(IVec& in,IVec& out,int items)
{
  int code = MPI_Alltoall(in.get_base(),items,MPI_INT,out.get_base(),items,
			  MPI_INT,MPI_COMM_WORLD);
  return code;
}

int MPIParMachine::all_to_all_displ(IVec&  in,IVec& sendcount,IVec& senddispl,
				     IVec& out,IVec& recvcount,IVec& recvdispl)
{
  int code = MPI_Alltoallv(in.get_base(),sendcount.get_base(),senddispl.get_base(),
			   MPI_INT, out.get_base(), recvcount.get_base(), 
			   recvdispl.get_base(), MPI_INT, MPI_COMM_WORLD);
  return code;
}

int MPIParMachine::all_to_all_displ(void*  in,IVec& sendcount,IVec& senddispl,
				     void* out,IVec& recvcount,IVec& recvdispl)
{
  int code = MPI_Alltoallv(in,sendcount.get_base(),senddispl.get_base(),
			   MPI_CHAR, out, recvcount.get_base(), 
			   recvdispl.get_base(), MPI_CHAR, MPI_COMM_WORLD);
  return code;
}

void MPIParMachine::barrier(const String& msg)
{
  if (my_id == 0)
    cout << " +++ Barrier: " << msg << " ...";
  cout.flush();
  MPI_Barrier(CompDomain);
  if (my_id == 0)
    cout << " reached " << endl;
  cout.flush();
}

int  MPIParMachine::nb_probe (int& flag)
{
  flag = 0;
  MPI_Status status;
  int code = MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,
			&flag,&status);
  source_id = status.MPI_SOURCE;
  msg_tag   = status.MPI_TAG;
  errcode   = status.MPI_ERROR;
  return code;
}

int  MPIParMachine::nb_send  (void* buf,int count,int dest,int tag)
{
  int code;
  int sum  = 0;
  int open = 0;
  MPI_Status status;
  int flag;
  for(int i = 0; i < rq.size(); i++)
    if (rq[i])
    {
      sum++;
      code = MPI_Test(rq[i],&flag,&status);
      if (!flag)
	open++;
    }
  
  cout << " Node: " << pm -> id() << " outstanding: " 
       << sum << " open: " << open <<  endl;
  
  if (rq[dest])
  {
    code = MPI_Test(rq[dest],&flag,&status);
    (ostream&) cout << code << " " << status.MPI_SOURCE << " " 
	 << status.MPI_ERROR << " " << flag << endl;
    if (!flag)
    {
      cout << " XXXX could not send " << endl;
      return NB_CANT_SEND;
    }
  }
  else
    rq[dest] = new MPI_Request;
    
  code = MPI_Issend(buf,count,MPI_CHAR,dest,tag,MPI_COMM_WORLD,
		   rq[dest]); 

  code = MPI_Test(rq[dest],&flag,&status);
  (ostream&) cout << code << " " << status.MPI_SOURCE << " " 
       << status.MPI_ERROR << " " << flag << endl;
 
  return NB_SENT;
}

int  MPIParMachine::receive  (void* buf,int& count)
{
  MPI_Status status;
  int code = MPI_Recv(buf,count,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,
		      MPI_COMM_WORLD,&status);
  MPI_Get_count(&status,MPI_CHAR,&count);
  source_id = status.MPI_SOURCE;
  msg_tag   = status.MPI_TAG;
  errcode   = status.MPI_ERROR;
  //cout << " RCV: " << id() << " code: " << code << " # " << count 
  // << " SRC: " << source_id 
  // << " TAG: " << msg_tag << " Err: " << errcode << endl;
  return code;
}

int  MPIParMachine::all_sum  (IVec& src,IVec& tar)
{
  assert(src.size() == tar.size());
  int code = MPI_Allreduce(src.get_base(),tar.get_base(),src.size(),
			   MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::all_sum  (int src,int& tar)
{
  int code = MPI_Allreduce(&src,&tar,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::all_sum  (double src,double& tar)
{
  int code = MPI_Allreduce(&src,&tar,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::all_sum  (Mat& src,Mat& tar)
{
  assert(src.size1() == tar.size1());
  assert(src.size2() == tar.size2());
  int code = MPI_Allreduce(src.get_base(),tar.get_base(),
			   src.size1()*src.size2(),
			   MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::all_sum  (NUMARRAY2D<float>& src,NUMARRAY2D<float>& tar)
{
  assert(src.size1() == tar.size1());
  assert(src.size2() == tar.size2());
  int code = MPI_Allreduce(src.get_base(),tar.get_base(),
			   src.size1()*src.size2(),
			   MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::broadcast(char* buf,int& count,int sender)
{
  int code;
  code = MPI_Bcast(&count,1,MPI_INT,sender,MPI_COMM_WORLD);
  code = MPI_Bcast(buf,count,MPI_CHAR,sender,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::broadcast(void* buf,int count,int sender)
{
  int code;
  code = MPI_Bcast(buf,count,MPI_CHAR,sender,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::broadcast(int& n,int sender)
{
  int code;
  code = MPI_Bcast(&n,1,MPI_INT,sender,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::broadcast(Mat& mat,int sender)
{
  int buffer[4];
  buffer[0] = mat.size1();
  buffer[1] = mat.size2();
  buffer[2] = mat.offset1();
  buffer[3] = mat.offset2();
  
  int code = MPI_Bcast(&buffer[0],4,MPI_INT,sender,MPI_COMM_WORLD);
  if (code != 0)
  {
    cout << " Node: " << id() << " Error in broadcast:: size of matrix " 
	 << code << endl; 
    error("Error in broadcast:: size of matrix");
  }
  
  if (id() != sender)
    mat.reset(buffer[0],buffer[1],buffer[2],buffer[3]);
  
  code = MPI_Bcast(mat.get_base(),mat.size1()*mat.size2(),
		   MPI_DOUBLE,sender,MPI_COMM_WORLD);
  return code;
}

int  MPIParMachine::broadcast(NUMARRAY2D<float>& mat,int sender)
{
  int buffer[4];
  buffer[0] = mat.size1();
  buffer[1] = mat.size2();
  buffer[2] = mat.offset1();
  buffer[3] = mat.offset2();
  
  int code = MPI_Bcast(&buffer[0],4,MPI_INT,sender,MPI_COMM_WORLD);
  if (code != 0)
  {
    cout << " Node: " << id() << " Error in broadcast:: size of matrix " 
	 << code << endl; 
    error("Error in broadcast:: size of matrix");
  }
  
  if (id() != sender)
    mat.reset(buffer[0],buffer[1],buffer[2],buffer[3]);
  
  code = MPI_Bcast(mat.get_base(),mat.size1()*mat.size2(),
		   MPI_FLOAT,sender,MPI_COMM_WORLD);
  return code;
}

int MPIParMachine::max_loc(double val)
{
  struct {
      double v;
      int    r;
  } in, out;
  in.r = pm->id();
  in.v = val;
  MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);  
  return out.r;
}
    
int MPIParMachine::min_loc(double val)
{
  struct {
      double v;
      int    r;
  } in, out;
  in.r = pm->id();
  in.v = val;
  MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);  
  return out.r;
}

int MPIParMachine::max_loc(int val)
{
  struct {
      int    v;
      int    r;
  } in, out;
  in.r = pm->id();
  in.v = val;
  MPI_Allreduce(&in,&out,1,MPI_2INT,MPI_MAXLOC,MPI_COMM_WORLD);  
  return out.r;
}

int MPIParMachine::min_loc(int val)
{
  struct {
      int    v;
      int    r;
  } in, out;
  in.r = pm->id();
  in.v = val;
  MPI_Allreduce(&in,&out,1,MPI_2INT,MPI_MINLOC,MPI_COMM_WORLD);  
  return out.r;
}

double MPIParMachine::max(double& val)
{
  double buf;
  MPI_Allreduce(&val,&buf,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
  return buf;
}
  
double MPIParMachine::min(double& val)
{
  double buf;
  MPI_Allreduce(&val,&buf,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
  return buf;
}

int MPIParMachine::max(int& val)
{
  int buf;
  MPI_Allreduce(&val,&buf,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);  
  return buf;
}

int MPIParMachine::min(int& val)
{
  int buf;
  MPI_Allreduce(&val,&buf,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);  
  return buf;
}

#endif




