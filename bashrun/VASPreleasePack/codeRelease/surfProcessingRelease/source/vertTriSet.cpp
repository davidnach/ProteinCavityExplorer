/****************************************************************************
 * SurfProcessing
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 *
 * This file is part of SurfProcessing
 * 
 * SurfProcessing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * SurfProcessing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SurfProcessing.  If not, see <http://www.gnu.org/licenses/>.  
 *
 * File: vertTriSet.cpp
 *       Special purpose lookup table for triangle generation
 *
 * Written by 
 *       Brian Y. Chen <chen@lehigh.edu>
 *
 * WWW URL: http://cse.lehigh.edu/~chen/
 * Email: chen@lehigh.edu
 * Documentation can be found here: 
 * http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
 *
 ***************************************************************************/

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

vertTriSet::vertTriSet()
{
	size = 0;
	edgeRoot = alloc_set(SP_MAP);

}

vertTriSet::~vertTriSet()
{

	int i;

	int size = size_set(edgeRoot);
	set_t tempSet;

	for(i = 0; i<size; i++){
		tempSet = (set_t) ith_map_set(edgeRoot, i);

		free_set(tempSet);
	}
	
	free_set(edgeRoot);
}

void vertTriSet::addEdge(int vert, int triangle)
{
	///insert both ways  first:
	set_t localSet = (set_t) mapsto_set(edgeRoot, vert);
	if(localSet == NIL){
		localSet = alloc_set(0);
		edgeRoot = put_set(edgeRoot, vert);
		localSet = put_set(localSet, triangle);
		edgeRoot = associate_set(edgeRoot, vert, (ptr_t) localSet);
	}
	else{
		localSet = put_set(localSet, triangle);
		edgeRoot = associate_set(edgeRoot, vert, (ptr_t) localSet);
	}
}

set_t vertTriSet::getNeighbors(int vert)
{
	int localSize = 0;
	int i = 0;

	set_t result = alloc_set(0);
	set_t localSet = (set_t) mapsto_set(edgeRoot, vert);

	if(localSet != NIL){
		localSize = size_set(localSet);
		for(i = 0; i<localSize; i++){
			result = put_set(result, localSet[i]);
		}
	}
	else{
		free_set(result);
		result = NULL;
	}

	return result;
}



