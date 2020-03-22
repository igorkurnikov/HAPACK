#include <stream.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>

const char NL   = '\n';
const char EOS  = '\0';
const char SP   = ' ';

FILE* depfile;
static int width  = 67;
static int indent = 15;

void error(char* msg,char* arg = NULL)
{
  if (arg == NULL)
    cerr << " +++ " << msg << NL;
  else 
    cerr << " +++ " << msg << SP << arg << NL;
  exit(3);
}

int parse(FILE* f,char* string)
{
  char *s = string;
  
  while(!feof(f) && *s != EOS)
  {
    char c = fgetc(f);
    // cout << "Parse: " << c << NL;
    if (c != *(s++)) break; 
  }
  
  return (*s == EOS);
  

}

int find(FILE* f,char *string)
{
  char c;
  
  while(!feof(f))
  {
    while(!feof(f) && (c = fgetc(f)) != string[0] && c != '"')
      ;
    if (c == '"') // skip strings
    {
      while(!feof(f) && (c = fgetc(f)) != '"')
	;
      continue;
    }
          
    if (c == string[0])
    {
      if (string[1] == EOS) return 1;
      if (parse(f,string+1)) return 1;
    }
  }
  return 0;
  
}


void scan(char* name,int& chars)
{
  FILE* f = fopen(name,"r");
  if (!f)
  {
    cerr << "+++ File not found: " << name << NL;
    return;
  }
    
  int found = find(f,"#include");
  while (found)
  {
    char c;
    while(!feof(f)) 
    {
      c = fgetc(f);
      if (!isspace(c)) break;
    }
    if (c == '"')
    {
      char buf[100];
      int cnt = 0;
      while(!feof(f) && (c = fgetc(f)) != '"' && (c != '\n'))
      {	
	buf[cnt++] = c;
        if (cnt >= 100) 
          error("Buffer Exhausted -- filename too long.");
      }
      if (c == '\n')
	error("Newline in String.");
      if (c != '"')
	error("failed to find a valid filename.");
      buf[cnt] = EOS;

      chars += strlen(buf) + 1;      
      
      if (chars > width - 2)
      {
          
	fprintf(depfile," \\\n");
        cout << " \\\n";      
        for(int i=0; i < indent;i++)
	{
	  fputc(SP,depfile);
	  cout << SP;
	}
	chars = indent;
      }
      fprintf(depfile," %s",buf);      
      cout << buf << SP;
      scan(buf,chars);      
    }      
    found = find(f,"#include");
  }

  fclose(f);  
}


main(int argc,char** argv)
{
  depfile = fopen("depends","w");
  
  argc--;

  while(argc > 0)
  {
    // strip suffix from the argument 
    char *name = argv[argc--];
    char *p = name;

    while(*p != EOS && *p != '.') p++;
    if(*p == '.')
      *p = EOS;

    char object[100];
    strcpy(object,name);
    strcat(object,".o");
    fprintf(depfile,"%s:",object);      
    cout << object << ":";
    for(int i=strlen(name)+2; i< indent;i++)
    {
      fputc(SP,depfile);
      cout << SP;
    }
    cout << name << ".c ";
    fprintf(depfile,"%s.c ",name);      
    int chrs = indent + strlen(name)+2;
    
    char filename[100];
    strcpy(filename,name);
    strcat(filename,".c");
    scan(filename,chrs);
    cout << NL;
    fputc(NL,depfile);
  }
  
  
  fclose(depfile);
  
    
    

}

