/*TEX
%
% IPACK - 
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
% with this package or write to: 
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: output.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Constructor}

The constructor takes a filename and an IO mode flag as arguments and
opens the file.  It also initializes and dumps -- or reads in and checks --
the first buffer.  Failures lead to abort. 
*/ 

#include "io_incl.h"
#include <errno.h>

#include "output.h"

OutputBuffer::OutputBuffer( const String& fname, const String& get_or_put ) {

  //  Open file in relevant mode.

  name     = fname;
  char *fs = fname();
  if ( get_or_put == "get" ) {
    f = fopen( fs, "r+b" );
    io_mode = 1;
  }
  if ( get_or_put == "put" ) {
    f = fopen( fs, "w+b" );
    io_mode = 0;
  }
  delete fs;

  //  Check file is opened
  
  if(!f) {
    cerr << "OutputBuffer::OutputBuffer() : Failure Opening IO-file: "
         << fname <<" "<< get_or_put;
    abort();
  }

  //  Read in or dump first buffer.

  if ( io_mode ) {

    read();
    if ( buf.compind[0] != 8 ) {
      cerr << "OutputBuffer::OutputBuffer() : buf.compind[0] != 8 : " 
           << buf.compind[0];
      abort();
    }

  }
  else {

    buf.no = 999;
    buf.compind[0] = 8;
    write();

  }
}

/*TEX
\subsection{Destructor}
The destructor closes the IO file after either checking it has
finished reading (this is not fatal, but a warning is issued), or
after writing a \name{twoel\_buffer} with $no = -1$ at the end of the file.
*/
OutputBuffer::~OutputBuffer() {

  if ( io_mode == 1 ) {
    if ( buf.no != -1 ) 
      cout << "Warning : OutputBuffer::~OutputBuffer() : buf.no != -1 : " 
           << buf.no << endl;
  }
  else if ( io_mode == 0 ) {
    if ( buf.no > 0 ) write();
    buf.no = -1;
    write();
  }
  fclose(f);

}

/*TEX
\subsection{Write and Read}
Write or read a buffer, and reset or return.
Yell on failure. 
*/
void OutputBuffer::write() {

  int code = fwrite( (char*) &buf, sizeof(buf), 1, f );
  if ( code != 1 ) {
    cerr << " ++++ OutputBuffer::write() : error on outputfile: " << name << endl;
    cerr.flush();
    abort();
  }
  buf.no = 0;

} 

int OutputBuffer::read() {

  int code = fread( (char*) &buf, sizeof(buf), 1, f );
  if ( code != 1 ) {
    cerr << " ++++ OutputBuffer::read() : error on inputfile: " << name << endl;
    cerr.flush();
    abort();
  }
  return buf.no;

} 
  
/*TEX
\subsection{Recieving information}
Gets and copies a buffer -- it might well be that deep copying is
not necessary in the long run, in which case change this to shallow
copying.
*/
/*
void recieve_data( const twoel_buffer* in_buf );
void OutputBuffer::recieve_data( const twoel_buffer* in_buf ) {

  buf.no = in_buf->no;
  int i;
  for ( i=0; i<buf.no; i++ ) {
    buf.value[i] = in_buf->value[i];
    buf.compind[i] = in_buf->compind[i];
  }

}
*/
