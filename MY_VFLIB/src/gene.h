/*----------------------------------------------------
 * gene.h
 * Interface of gene.cc
 * Random generation of isomorphic ARGraphs
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: gene.h,v 1.1.1.1 2008/04/08 19:44:27 igor Exp $
 ----------------------------------------------------*/

/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: gene.h,v $
 *   Revision 1.1.1.1  2008/04/08 19:44:27  igor
 *
 *
 *   Revision 1.1.1.1  2003/10/15 06:24:11  igor
 *   load vflib
 *
 *   Revision 1.2  1998/12/08 13:31:01  foggia
 *   Minor changes
 *
 *---------------------------------------------------------------*/

#ifndef GENE_H

#include "argraph.h"



void Generate(int nodes, int edges, Graph **g1, Graph **g2, 
              bool connected=true);

#endif
