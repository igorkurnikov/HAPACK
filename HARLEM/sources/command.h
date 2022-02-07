/*! \file command.h

   Classes to process RASMOL-like and PYTHON text command  
   
   \author Igor Kurnikov
   \date 1997-2003
*/

#if !defined(COMMAND_H)
#define COMMAND_H

#include "haio.h"
#include "haconst.h"
#include "hastl.h"
#include "hastring.h"

/* Format values are related to Tokens */
#define Tok2Format(x) ((x)-359)
#define Format2Tok(x) ((x)+359)

const int FormatPDB      = 1;
const int FormatMacroMod = 2;
const int FormatGaussian = 3;
const int FormatAlchemy  = 4;
const int FormatNMRPDB   = 5;
const int FormatCharmm   = 6;
const int FormatBiosym   = 7;
const int FormatMOPAC    = 8;
const int FormatSHELX    = 9;
const int FormatMol2     = 10;
const int FormatFDAT     = 11;
const int FormatMMDB     = 12;
const int FormatMDL      = 13;
const int FormatXYZ      = 14;
const int FormatCIF      = 15;
const int FormatCEX      = 16;
const int FormatRWF      = 17;
const int FormatHarlem    = 18;
const int FormatAmberPrep = 19;
const int FormatAmberTop = 20;
const int FormatAmberOff = 21;
const int FormatHIN      = 22;
const int FormatNRG      = 23;

const int FormatGUESS    = 0;

class AtomExpr;
class MolSet;

class CmdParser
//! Class to process RASMOL-like text commands
{
public:
	CmdParser();
	virtual ~CmdParser();
	
	int CurToken;                //!< the integer value of the token in Keywords map

	deque<std::string> cmd_history;  //!< History of commands entered
	static int max_history_save;
	int cur_history_cmd;
	
	static StrIntMap Keywords;        //!< Map of Keywords
	static int InitKeywords();                                 //!< function to initiate map of keywords
	static int RegisterKeyword(std::string keyw, const int itok); //!< Register(add) keyword to the Keywords map 

	int TokenValue;          //!< Integer value of the current token
	double TokenValueFloat;  //!< Float value of the current token
	std::string TokenIdent;  //!< String value of the current token (if in "" quotes)
	size_t str_start_pos;       //!< Starting position of the string token (if without "" quotes)
	size_t cursor_pos;       //!< Current cursor position in the command line 
	                                                        
	int SetCmdLine(const std::string& cmd_line); //!< Set Command Line
	const char* GetCmdLine();  //!< Get current Command Line
	const char* GetStartPosSubstr(); //!< Get Substring of the command line staring with str_start_pos
	std::string RollHistory(int step);    //!< Retrieve a command from the history

	void ResetCursorPosition(); //!< Set Cursor position at the beginning of the command line
	void CommandError(const char* error );                        
	int LookUpKeyword();                     //!< Get an integer value for a current token(TokenIdent) return IdentTok if keyword is not registered in Keywords map
	int FetchToken();                        //!< Get Next token in the line after CursorPtr
	int NextIf(int token, const char* error );
	int ParseColour(int& RVal, int& GVal, int& BVal); //!< Set RGB value corresponging to the colour
	AtomExpr* ParseRange(int neg );
	AtomExpr* ParseExpression(int level, MolSet* pmset); //!< Parse str_parse and form a logical expression in QParse

protected:
	std::string CurLine;    //!< Line with the command to process

};


#ifdef _WIN32
#define DirChar  '\\'
#else
#define DirChar  '/'
#endif


#endif /* !defined(COMMAND_H) */
