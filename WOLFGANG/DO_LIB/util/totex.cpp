#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

const char EOS = '\0';
const char  NL = '\n';
const char  SP = ' ';

void error(char* msg)
{
  cerr << " +++ Fatal Error: " << msg << NL;
  abort();
}

const int BUFLEN   = 1024;
const int NOINFILE = 1;
const int SUCCESS  = 0;  

class ToTeX
{

  char* in_fname;
  char* out_fname;
  char* outbuf;
  char* inbuf;
  char  delim;
  
  int found_tab;  // delimiter of last token was a tab !
  
  ifstream infile;
  ofstream outfile;

  void  strip_suffix(char* name);
  char* process_token(char* intoken);
  int   read_token();  

  enum  TEX_MODES { OFF, ON };
    
  int   tex_mode;
  void  texmode_on();
  void  texmode_off();
  
public:


    ToTeX(char* infilename);
   ~ToTeX();
  
int process(char* name); 
};

ToTeX::ToTeX(char* infilename)
{
  found_tab = 0;
  
  outbuf    = new char[BUFLEN];
  inbuf     = new char[BUFLEN];
  
  in_fname  = new char[strlen(infilename) + 10];
  strcpy(in_fname,infilename);
  strip_suffix(in_fname);

  out_fname = new char[strlen(in_fname) + 10];
  strcpy(out_fname,in_fname);
  strcat(out_fname,".tex");    
    
  cout << " +++ Opening Output File: " << out_fname << NL;
  outfile.open(out_fname,ios::out);
  if (!outfile)
    error("Could not open output file.");

  // attempt to process the .h file

  int len = strlen(in_fname);
  in_fname[len  ] = '.';
  in_fname[len+1] = 'h';
  in_fname[len+2] = EOS;

  tex_mode = ON;
  process(in_fname);   

  // attempt to process the .h file
  in_fname[len+1] = 'c';

  tex_mode = ON;
  process(in_fname);   

}

ToTeX::~ToTeX()
{
  delete in_fname;
  delete out_fname;
  delete outbuf;  
}

char* ToTeX::process_token(char* intoken)
{
  if (strcmp(intoken,"/*TEX")==0) 
  {
    texmode_on();
    intoken[0] = EOS;
  }

  if (strcmp(intoken,"\t")==0) 
  {
    strcpy(intoken,"        ");
    return intoken;
  }
  
  if (strcmp(intoken,"*/")==0 && tex_mode) 
  {
    texmode_off();
    intoken[0] = EOS;
  }
  return intoken;  
}

void ToTeX::texmode_on()
{
  if (tex_mode == ON)
    error(" TeXmode is already on!");
  outfile << NL << "\\end{verbatim} " << NL
          << "{\\ \\hrulefill End of Code \\hrulefill} \\ " 
    << "\\normalsize" << endl;
  tex_mode = ON;
}

void ToTeX::texmode_off()
{
  if (tex_mode == OFF)
    error(" TeXmode is already off!");
  tex_mode = OFF;
  outfile << NL << "\\vspace{0.1em} \\small " << endl 
          << "{\\ \\hrulefill Start of Code \\hrulefill} " 
	  << endl << "\\begin{verbatim} " << NL;
}

int  ToTeX::read_token()
{
  if (found_tab) // if last token was tab ....
  {
    found_tab = 0;
    strcpy(inbuf,"\t");
    return 1;
  }
  
  int count = 0;
  char c;
  
  while(infile.get(c))
  {
    // check for delimiters
    delim = c;
    if (isspace(c)) 
    {
      if (c=='\t') found_tab = 1; // set the flag that tab was found
      inbuf[count] = EOS;
      return 1;
    }
    
    // add char to buffer
  
    inbuf[count++] = c;    
                   
    if (count == BUFLEN-1)
    {
      cerr << " +++ Warning: Buffer Overflow." << NL;
      delim = EOS;
      inbuf[count] = EOS;
      return 1;
    }
  }
  return 0;
}


int  ToTeX::process(char* name)
{
  cout << " +++ Opening Input File: " << name << NL;
  infile.open(name,ios::in);
  if (!infile) return NOINFILE;

  outfile << "\\newpage " << NL;
  outfile << "\\markboth{" << name << "}{" << name << "}" <<  NL;

  if (read_token())
  {
    if (strcmp(inbuf,"/*TEX")==0)
      tex_mode = ON;
    else
    {
      if (tex_mode == ON)
	texmode_off();
      outfile << process_token(inbuf);  
    }
  }
  
  while(read_token())
  {
    outfile << process_token(inbuf);  
    if (delim != EOS)
      outfile << delim;
  }
  
  infile.close();  

  if (tex_mode == OFF)
    outfile << "\\end{verbatim} " << endl; 

  return SUCCESS;
  
}


void ToTeX::strip_suffix(char* name)
{
  cout << "NAME: " << name << endl;
  
  for(int i = strlen(name)-1; i>0; i--)
    if (name[i] == '.')
    {
      name[i] = EOS;
      break;
    }
  
  cout << "stripped name: " << name << NL;  
}


main(int argc,char** argv)
{
  if (argc < 2)
    error("Usage: totex <filename>");
  
  ToTeX convert(argv[1]);
  return 0;  
}
