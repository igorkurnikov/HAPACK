/*------------------------------------------------------------------
 * match.h
 * Header of match.cc
 * Declaration of the match function
 *
 * Author: P. Foggia
 * $Id: match.h,v 1.1.1.1 2008/04/08 19:44:27 igor Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: match.h,v $
 *   Revision 1.1.1.1  2008/04/08 19:44:27  igor
 *
 *
 *   Revision 1.1.1.1  2003/10/15 06:24:11  igor
 *   load vflib
 *
 *   Revision 1.1  1998/09/29 09:49:48  foggia
 *   Initial revision
 *
 *----------------------------------------------------------------*/


#ifndef MATCH_H
#define MATCH_H

#include "argraph.h"
#include "state.h"

/*------------------------------------------------------------
 * Definition of the match_visitor type
 * a match visitor is a function that is invoked for
 * each match that has been found.
 * If the function returns false, then the next match is
 * searched; else the seach process terminates.
 -----------------------------------------------------------*/
typedef bool (*match_visitor)(int n, node_id c1[], node_id c2[], 
                              void *usr_data);

bool match(State *s0, int *pn, node_id c1[], node_id c2[]);

int match(State *s0, match_visitor vis, void *usr_data=NULL);

#endif
