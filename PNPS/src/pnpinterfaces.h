//
// C++ Interface: pnpinterfaces
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PNPINTERFACES_H
#define PNPINTERFACES_H
class Data;
/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class GenericSolver{
  public:
    GenericSolver();

    virtual ~GenericSolver();
    virtual int InitZero()=0;//!<initialize zero values
    virtual int Clear()=0;//!< delete all internal arrays and objects
    virtual int ShowParameters()=0;//!<show parameters of the solver
    //!show properies can be diffrent from parameters of the solver, valid after InitSolver()
    virtual int ShowProperties()=0;
    
    virtual int InitSolver()=0;//!<initiate internal arrays
    virtual int Solve()=0;//!< Solve problem

};

#endif
