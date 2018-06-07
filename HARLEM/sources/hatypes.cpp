/*! \file hatypes.cpp

    Some basic types to use in HARLEM 
 
    \author Igor Kurnikov  
    \date 2010-
*/

#include <mpi.h>

#include "hatypes.h"

int HaEnum::Bcast(MPI_Comm &comm, int root)
{
	return MPI_Bcast(&this->value(),1,MPI_INT,root,comm);	
}

std::vector<int> HaEnum1::all_values() const
{
	IntStrMap& lbl_map = GetLabelsMap();
	std::vector<int> values;
	IntStrMap::iterator itr;
	for(itr = lbl_map.begin(); itr != lbl_map.end(); itr++)
	{
		values.push_back((*itr).first);
	}
	return values;
}


std::vector<std::string> HaEnum1::GetAllLabels() const
{
	IntStrMap& lbl_map = GetLabelsMap();
	std::vector<std::string> labels;
	IntStrMap::iterator itr;
	for(itr = lbl_map.begin(); itr != lbl_map.end(); itr++)
	{
		labels.push_back((*itr).second);
	}
	return labels;
}
