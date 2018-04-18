/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 *
 * This file is part of VASP
 * 
 * VASP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * VASP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with VASP.  If not, see <http://www.gnu.org/licenses/>.  
 *
 * File: Lattice.cpp
 *       A data structure for managing marching cubes
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

// #####################################################################################################################################
// #####################################################################################################################################
// ####################    ###########    ######            ##            ##        #####       #####         ##########################
// ####################    ##########      #####            ##            ##        ###           ###         ##########################
// ####################    #########        ########    ##########    ########    ####    #####    ##    ###############################
// ####################    ########    ##    #######    ##########    ########    ####    ###########    ###############################
// ####################    #######    ####    ######    ##########    ########    ####    ###########       ############################
// ####################    #######    ####    ######    ##########    ########    ####    ###########       ############################
// ####################    #######            ######    ##########    ########    ####    ###########    ###############################
// ####################    #######            ######    ##########    ########    ####    ###########    ###############################
// ####################    #######    ####    ######    ##########    ########    ####    #####    ##    ###############################
// ####################         ##    ####    ######    ##########    ######        ###           ###         ##########################
// ####################         ##    ####    ######    ##########    ######        #####       #####         ##########################
// #####################################################################################################################################
// #####################################################################################################################################

// this code uses the following cube layout, with 3D axes to the right.
//       3---7            ^
//       |\  |\       \   |
//       | 1---5       +y |
//       | | | |        \ z+
//       2-|-6 |         \|  
//        \|  \|          0--x+-->    
//         0---4

#include "Lattice.h"

/////////////////////////////////////////////////////////////////////////////////////
// Nothing is computed in the constructor.  Just set up the data.
/////////////////////////////////////////////////////////////////////////////////////
Lattice::Lattice( SurfaceObject * s )
{
	
	surf = s;							// surface data
	res = 0;							// null resolution
	
	///set the initial boundary values
	xneg = HUGE_VAL;	xpos = -HUGE_VAL;
	yneg = HUGE_VAL;	ypos = -HUGE_VAL;
	zneg = HUGE_VAL;	zpos = -HUGE_VAL;

	///set the initial grid extents.
	xdim = 0;
	ydim = 0;
	zdim = 0;

	nonempty = alloc_set(SP_MAP);			// Set pointing to cubes that contain geometry.
									// This array is indexed as described in the	
									// insertTriangles() comment above.

	//nonemptyGroups = alloc_set(0);		// A set of sets listing indices of connected 
									// components of adjacent non-empty cubes.
									// Allocated in classifyCubes(), so commented here

	//emptyGroups = alloc_set(0);			// A set of sets listing indices of connected 
									// components of adjacent empty cubes.
									// Allocated in classifyCubes(), so commented here

	neighborMatrix = alloc_set(SP_MAP);	// A set of sets of ints.  Top level indices correspond 
									// to EMPTY connected components and link to empty sets, 
									// whose (second-level) indices correspond to other empty
									// sets, and link to an int[1], which stores the index
									// of the cube that is adjacent to cubes of both empty
									// sets.
									
	edges = alloc_set(SP_MAP);			// This set contains sets that map to a point in space.
	edgePoints = alloc_set(SP_MAP);		// The top level set is indexed by points, and the
									// lower level set is also indexed by points, so you get
									// a way to index point-to-point connections - i.e. edges.
									// The lower level set contains double* values that are the
									// actual double positions, so that they are never recomputed.
									// - note that only one instance is stored, so if you
									// look for edges, you must enumerate the point with lower 
									// index, then the point with upper index.

	points = alloc_set(SP_MAP);			// This set contains pointers to point_inside or
	point_inside = new int[1];			// point_outside, indexed on point indices.  this simply
	point_inside[0] = LATTICE_INTERIOR;	// states if a given point is inside or outside, and is only
	point_outside = new int[1];			// applied to points that are not already set by cube	
	point_outside[0] = LATTICE_EXTERIOR;	// flagging.

	numGeometricChecks = 0;				// diagnostic info recording

}

/////////////////////////////////////////////////////////////////////////////////////
// destructor
/////////////////////////////////////////////////////////////////////////////////////
Lattice::~Lattice()
{
	int i = 0;
	int j = 0;

	// This is the surf file we use for initial data.  Do Not Delete it.
	//Do not delete because it is only referenced, not copied.
	//delete(surf);

	// This is a set that lists the cubes that have triangles.
	// So delete each set, then delete the whole set
	for(i = 0; i<size_set(nonempty); i++){
		set_t tempSet = (set_t) mapsto_set(nonempty, nonempty[i]);
		free_set(tempSet);
	}
	free_set(nonempty);

	// this is a set of sets containing the indices of adjacent cubes that contain triangles.
	// So delete each set, then delete the whole set
	for(i = 0; i<size_set(nonemptyGroups); i++){
		set_t tempSet = (set_t) mapsto_set(nonemptyGroups, nonemptyGroups[i]);
		free_set(tempSet);
	}
	free_set(nonemptyGroups);

	// this is a set of sets containing the indices of adjacent cubes that do not contain triangles.
	// So delete each set, then delete the whole set
	for(i = 0; i<size_set(emptyGroups); i++){
		set_t tempSet = (set_t) mapsto_set(emptyGroups, emptyGroups[i]);
		free_set(tempSet);
	}
	free_set(emptyGroups);

	// A set of sets of ints.  Delete carefully.
	// So delete each int, then delete the set that used to contain ints, then delete the whole set
	for(i = 0; i<size_set(neighborMatrix); i++){
		set_t tempSet = (set_t) mapsto_set(neighborMatrix, neighborMatrix[i]);
		for(j = 0; j<size_set(tempSet); j++){
			int * temp = (int *) mapsto_set(tempSet, tempSet[j]);
			delete[](temp);
		}
		free_set(tempSet);
	}
	free_set(neighborMatrix);

	//An array that stores the interior/exterior status of all connected component groups
	delete[](flags);

	// This is a set of sets indicating the connection between cornerIndexes and other corner indexes
	// So delete each int, then delete the set that used to contain ints, then delete the whole set
	for(i = 0; i<size_set(edges); i++){
		set_t thisSet = (set_t) mapsto_set(edges, edges[i]);
			for(j = 0; j<size_set(thisSet); j++){
				int * pt = (int *) mapsto_set(thisSet, thisSet[j]);
				delete[](pt);
			}
		free_set(thisSet);
	}
	free_set(edges);

	// This is a set of double*'s that are indexed by edges
	// So delete each double*, then delete the whole set
	for(i = 0; i<size_set(edgePoints); i++){
		double * pt = (double *) mapsto_set(edgePoints, edgePoints[i]);
		delete[](pt);
	}
	free_set(edgePoints);

	//this is a set that contains arrays of size one, but we dont need to delete them individually
	free_set(points);			///note that all the pointers inside point to one 
	delete[](point_inside);		///point_inside or point_outside, so jsut delete them
	delete[](point_outside);		///separately.
}



// ##################################################################################################################
// ##     ##     ##    #  #####  #### ####     ##  ####  ######     ###     ##  ###  #      #     ###     ##  #######
// ##  ##  #  ##  ##  ##   ###   ###   ###  ##  #  ####  #####  ###  #  ###  #   ##  ###  ###  ##  #  ###  #  #######
// ##  ##  #  ##  ##  ##    #    ##     ##  ##  #  ####  #####  ###  #  ###  #    #  ###  ###  ##  #  ###  #  #######
// ##  ##  #  ##  ##  ##  #   #  #   #   #  ##  ##  ##  ######  ######  ###  #  #    ###  ###  ##  #  ###  #  #######
// ##     ##    ####  ##  ## ##  #  ###  #    #####    #######  ######  ###  #  ##   ###  ###    ###  ###  #  #######
// ##  #####  #  ###  ##  #####  #       #  #  #####  ########  ######  ###  #  ###  ###  ###  #  ##  ###  #  #######
// ##  #####  ##  ##  ##  #####  #  ###  #  ##  ####  ########  ###  #  ###  #  ###  ###  ###  ##  #  ###  #  #######
// ##  #####  ##  ##  ##  #####  #  ###  #  ##  ####  ########  ###  #  ###  #  ###  ###  ###  ##  #  ###  #  #######
// ##  #####  ##  #    #  #####  #  ###  #  ##  ####  #########     ###     ##  ###  ###  ###  ##  ##     ##     ####
// ##################################################################################################################


/////////////////////////////////////////////////////////////////////////////////////
// setBounds() sets the boundaries of the Lattice, given the surface (surf)
// stored in the Lattice ONLY.  For CSG operations with multiple lattices,
// Lattices are synchronized for CSG using the synchronize() function.
/////////////////////////////////////////////////////////////////////////////////////
void Lattice::setBounds(double r)
{
	int i = 0;
	
	////////////////////////////////////////
	printf("STATE [setBounds()          ]: Exploring Surface Bounds... \n");
	////////////////////////////////////////

	// set the resolution
	res = r;

	///get Surface Bounds
	for(i = 0; i<surf->numPoints; i++){
		stdoutProgressBar(i, surf->numPoints);
		if(surf->surfacePoints[3*i+0] < xneg){ xneg = surf->surfacePoints[3*i+0]; }
		if(surf->surfacePoints[3*i+0] > xpos){ xpos = surf->surfacePoints[3*i+0]; }
		if(surf->surfacePoints[3*i+1] < yneg){ yneg = surf->surfacePoints[3*i+1]; }
		if(surf->surfacePoints[3*i+1] > ypos){ ypos = surf->surfacePoints[3*i+1]; }
		if(surf->surfacePoints[3*i+2] < zneg){ zneg = surf->surfacePoints[3*i+2]; }
		if(surf->surfacePoints[3*i+2] > zpos){ zpos = surf->surfacePoints[3*i+2]; }
	}

	////////////////////////////////////////
	printf("INFO  [Input Bounds         ]: xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	////////////////////////////////////////

	///random nudge to offset axis aligned geometry, + 1 row of low-side padding.
	double xnudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double ynudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double znudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	xneg -= (xnudge+res);	// We nudge to make the cube bigger.
	yneg -= (ynudge+res);	// This way there is no chance something does not fit.
	zneg -= (znudge+res);

	///+1 for fractional cube round up (int conversion truncates)
	///+1 for high side padding
	///= +2 total additions.
	xdim = ((int) ((xpos-xneg)/res)) + 2;       ///this is the total number of CUBES.  Not lines.
	ydim = ((int) ((ypos-yneg)/res)) + 2;
	zdim = ((int) ((zpos-zneg)/res)) + 2;

	///Reset the high dimensions based on the sizes of the cubes, since there is fractional roundup.
	xpos = xneg + xdim*res;
	ypos = yneg + ydim*res;
	zpos = zneg + zdim*res;
        
	////////////////////////////////////////
	printf("INFO  [Expanded Dimensions  ]: xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	printf("INFO  [Cube Dimensions      ]: xdim: %i, ydim: %i, zdim: %i\n", xdim, ydim, zdim);
	printf("INFO  [Lattice Resolution   ]: %f\n", res);
	printf("INFO  [Total Num. of Cubes  ]: %i\n", xdim*ydim*zdim);
	////////////////////////////////////////

}



/////////////////////////////////////////////////////////////////////////////////////
// Scan conversion: 
// This function inserts all triangles from the surfaceObject into cubes.
// The cubes are virtual.  They are conceptually indexed in a grid, with indices
//    i*(ydim*zdim) + j*(zdim) + k;
/////////////////////////////////////////////////////////////////////////////////////
void Lattice::insertTriangles()
{
	////////////////////////////////////////
	printf("INFO  [insertTriangles()    ]: Now Inserting %i Triangles from input\n", surf->numTriangles);
	printf("STATE [insertTriangles()    ]: Inserting Triangles into Lattice... ");
	////////////////////////////////////////
	int l = 0;
	for(l = 0; l<surf->numTriangles; l++){
		stdoutProgressBar(l, surf->numTriangles);
		int i = 0;	int j = 0;	int k = 0;
		
		// get the triangle we will work on.
		double * t = surf->getTriangle(l); 

		// Find the bounding box of the triangle.
		double sxmax = -HUGE_VAL;	double sxmin = HUGE_VAL;
		double symax = -HUGE_VAL;	double symin = HUGE_VAL;
		double szmax = -HUGE_VAL;	double szmin = HUGE_VAL;
		for(i = 0; i<3; i++){
			if(sxmax < t[3*i+0]){ sxmax = t[3*i+0]; }	if(sxmin > t[3*i+0]){ sxmin = t[3*i+0]; }
			if(symax < t[3*i+1]){ symax = t[3*i+1]; }	if(symin > t[3*i+1]){ symin = t[3*i+1]; }
			if(szmax < t[3*i+2]){ szmax = t[3*i+2]; }	if(szmin > t[3*i+2]){ szmin = t[3*i+2]; }
		}

		// map the indices as if the grid was infinite.
		int ihi, ilo, jhi, jlo, khi, klo;
		if( fmod( (sxmin-xneg), res ) == 0 ){ ilo = ((int)((sxmin-xneg)/res)) -1; }
		else{ ilo = ((int)((sxmin-xneg)/res)); }
		if( fmod( (symin-yneg), res ) == 0 ){ jlo = ((int)((symin-yneg)/res)) -1; }
		else{ jlo = ((int)((symin-yneg)/res)); }
		if( fmod( (szmin-zneg), res ) == 0 ){ klo = ((int)((szmin-zneg)/res)) -1; }
		else{ klo = ((int)((szmin-zneg)/res)); }
		ihi = (int)( (sxmax-xneg)/res );	// if it is on the exact high end of a cube, it will 
		jhi = (int)( (symax-yneg)/res );	//  give the whole additional cube, so the checks above
		khi = (int)( (szmax-zneg)/res );	//  do not apply to ihi, jhi, khi.

		// Check index bounds.  Breaking index bounds should NOT be possible
		// Unless the lattice was not built for this surface.
		if( ihi>=xdim || ilo<0 || jhi>=ydim || jlo<0 || khi>=zdim || klo<0 ){
			printf("ERROR [Lattice Bounds]: Attempted triangle insertion violates Lattice bounds.\n");
			printf("ERROR [Lattice Bounds]: This Lattice was not built using your surface.  Exitting.\n");
			exit(1);
		}

		// retrieve the cubes.  Note that the cubes are indexed in three dimensions, so we must
		// first span the indices, then transform them into a single dimension index (see myindex)
		for(i = ilo; i<=ihi; i++){
			for(j = jlo; j<=jhi; j++){
				for(k = klo; k<=khi; k++){
					// get the index
					int myIndex = (i*ydim*zdim) + (j*zdim) + k;
					
					// make sure the triangle actually intersects the cube.
					double xn = xneg + ((i+0)*res);
					double xp = xneg + ((i+1)*res);
					double yn = yneg + ((j+0)*res);
					double yp = yneg + ((j+1)*res);
					double zn = zneg + ((k+0)*res);
					double zp = zneg + ((k+1)*res);

					// This function decides if a given triangle intersects a given cube.
					// If it doesnt, move on.  If it does, store it's index.					
					// FIXME: this function call could possibly be made faster.
					if( !triCubeIntersect(t, xn, xp, yn, yp, zn, zp) ){ 
						continue; 
					}
					
					// find out of this cube has been made yet
					// We do not allocate any cubes unless somethign must be stored in them
					set_t mySet;

					// if it was made, get it
					if( contains_set(nonempty, myIndex) ){ mySet = (set_t) mapsto_set(nonempty, myIndex); }

					// if it was not made, allocate it.
					else{ mySet = alloc_set(0); }
					
					// store the index of the triangle in the cube.
					mySet = put_set( mySet, l );

				//	if(index1 == 4026 || index2 == 4289 || index3 == 4027){
				//		printf("triangle %i, [%i %i %i] inserted into cube %i\n", l, index1, index2, index3, myIndex);
				//	}

					// store the cube into the nonempty cubes, since we have added to the
					// set and its pointer may have changed on set expansion.
					nonempty = associate_set( nonempty, myIndex, mySet );
				}
			}
		}

		//clear the triangle, since we just need the index
		delete[](t);
	}

//	////////////////////////////////////////
//	printf("INFO  [insertTriangles()]: Done Inserting Triangles... \n");
//	////////////////////////////////////////
}



/////////////////////////////////////////////////////////////////////////////////////
// Cube sign generation step 1
// This function finds all the connected components of empty and nonempty cubes.
// Lists of indices of empty cubes go into emptyGroups
// Lists of indices of nonempty cubes go into nonemptyGroups
// Allocates the flags array, which could not be allocated without the number of
// empty groups
/////////////////////////////////////////////////////////////////////////////////////
void Lattice::classifyCubes()
{
	int i = 0;
	set_t empty = alloc_set(0);

	///fill a temporary array with cubes that are empty.
	for(i = 0; i<xdim*ydim*zdim; i++){
		if( !contains_set(nonempty, i) ){
			empty = put_set(empty, i);
		}
	}
	
	///these are the classes of cubes.
	////////////////////////////////////////
	printf("STATE [classifyCubes()      ]: Finding connected non-empty cubes... ");
	////////////////////////////////////////
	nonemptyGroups = unionFind( empty );	//THis is a list of lists of adajcent **NONEMPTY** cubes
	////////////////////////////////////////
	printf("INFO  [# Nonempty Groups    ]: %i\n", size_set(nonemptyGroups) );
	////////////////////////////////////////

	////////////////////////////////////////
	printf("STATE [classifyCubes()      ]: Finding connected empty cubes... ");
	////////////////////////////////////////
	emptyGroups = unionFind( nonempty );	//THis is a list of lists of adajcent **EMPTY** cubes
	////////////////////////////////////////
	printf("INFO  [# Empty Cube Groups  ]: %i\n", size_set(emptyGroups) );
	////////////////////////////////////////

	free_set(empty);
	
	// allocate the flags array. These flags will record the interior/exterior
	// status of groups of adjacent empty cubes.
	flags = new int[size_set(emptyGroups)];
	for(i = 0; i<size_set(emptyGroups); i++){
		flags[i] = LATTICE_NOT_SET;
	}
}




/////////////////////////////////////////////////////////////////////////////////////
// Cube sign generation step 2
// This function fills the neighborMatrix, with segment-based testing on
// the **nonempty** cubesets.  In other words, we are finding groups of empty cubes
// that are adjacent to each other - i.e. on the other side of **nonempty** cubes
// from each other.
/////////////////////////////////////////////////////////////////////////////////////
void Lattice::fillNeighboringSets()
{
	int i = 0;
	int j = 0;
	int k = 0;

	////////////////////////////////////////
	printf("STATE [fillNeighboringSets()]: Discovering lattice topology via filled-cube adjacency... ");
	////////////////////////////////////////
	//we are going through every cube here, looking for neighbors that have been classified already.
	for(i = 0; i<size_set(nonempty); i++){
		stdoutProgressBar(i, size_set(nonempty));
		// get the cube
		int thisNEcube = nonempty[i];
		// find neighboring empty groups
		set_t neighborGroups = getNeighboringEmptyGroups(thisNEcube);
		
		// generate pairs to put into neighborMatrix
		// note that if the cube neighbors only one group, this automatically has 
		// no effect, as desired.  If the cube has more than 2, this works fine also.
		//////
		//	printf("This many neighboring Groups: %i\n", size_set(neighborGroups) );
		// Here there may only be one kind of neighboring group, because you are on a tiny
		// non-empty spot that is not large enough to surround an empty space, only have
		// some non-empty space alone somewhere.
		// **THUS - it can be the case that nothing is ever entered into neighborMatrix**
		for(j = 0; j<size_set(neighborGroups); j++){
			// the firstIndex is the first half of the pair
			int firstIndex = neighborGroups[j];
			for(k = j+1; k<size_set(neighborGroups); k++){
				// the secondIndex is the second half of the pair
				int secondIndex = neighborGroups[k];
				
				// now insert forwards (i.e. first-second).
				set_t subSet;
				if( contains_set(neighborMatrix, firstIndex) ){
					subSet = (set_t) mapsto_set(neighborMatrix, firstIndex);
				}
				else{ subSet = alloc_set(SP_MAP); }
				// if it already contains first and second, then we do not need to insert.
				// obviuosly if we had to allocate the subSet, then it does not contain it though.
				if( contains_set(subSet, secondIndex) ){ break; }
				else{
					// Now store the cube where we found the adjacent groups.
					int * box = new int[1];
					box[0] = thisNEcube;
					subSet = associate_set(subSet, secondIndex, box);
					neighborMatrix = associate_set(neighborMatrix, firstIndex, subSet);
				}

				// now insert backwards (i.e. second-first).
				if( contains_set(neighborMatrix, secondIndex) ){
					subSet = (set_t) mapsto_set(neighborMatrix, secondIndex);
				}
				else{ subSet = alloc_set(SP_MAP); }
				// if it already contains second and first, then we do not need to insert.
				if( contains_set(subSet, firstIndex) ){ break; }
				else{
					//now store the cube where we found the adjacent gorups.
					int * box = new int[1];
					box[0] = thisNEcube;
					subSet = associate_set(subSet, firstIndex, box);
					neighborMatrix = associate_set(neighborMatrix, secondIndex, subSet);
				}
			}
		}
		//free this now that we are done with it.
		free_set(neighborGroups);
	}

	////////////////////////////////////////
//	printf("STATE [fillNeighboringSets()]: Lattice topology exploration complete.\n");
	////////////////////////////////////////
}







/////////////////////////////////////////////////////////////////////////////////////
// Cube sign generation step 3 
// This function runs segment-based tests on adjacent empty groups, based on the information
// in the neighborMatrix.  The neighbormatrix specifies specific cubes that are adjacent to
// two different groups of empty connected cubes.  Segment tests are run in the vicinity of
// those specific cubes, to ensure that the segments tested are as short as possible.
// Furthermore, using breath first search, tests are only run as needed to propagate inside-
// outside information from the exterior (which is, by definition, exterior) too all other
// empty groups.  As a result except for some special cases, all empty groups should be
// properly flagged in using this function.  flagProximalSets covers all remaining cases.
/////////////////////////////////////////////////////////////////////////////////////
void Lattice::populateFlags()
{
	////////////////////////////////////////
	printf("STATE [populateFlags()      ]: Assigning interior/exterior state via neighbor segment testing...\n");
	////////////////////////////////////////
	int i = 0;
	int j = 0;

	// first, find the cubeSet that contains 0, 
	// the cube that is empty by definition (because it is the exterior)
	int containsZero = -1;
	for(i = 0; i<size_set(emptyGroups); i++){
		set_t thisSet = (set_t) mapsto_set(emptyGroups, i);
		if( contains_set(thisSet, 0) ){
			containsZero = i;
			break;
		}
	}
	flags[containsZero] = LATTICE_EXTERIOR;	//this is true by definition.
	
	//error check.
	if( containsZero == -1 ){
		printf("ERROR [populateFlags()]: No empty set contains the cube index 0.  Contradiction!\n");
		exit(1);
	}
	
	//simplicity check.  If there is only one exterior, then we have already set everything,
	//and the neighborMatrix contains no neighbors, so if we continue, we will hit a null pointer and seg fault.
	if( size_set(emptyGroups) == 1 ){
		return;
	}

	//allocate the working set for our breadth first search
	set_t workingSet = alloc_set(0);
	workingSet = put_set(workingSet, containsZero);

	//now, as long as the workingSet is nonEmpty, we continue our breadth first search.
	while( size_set(workingSet)>0 ){
		//nextSet is the set that will contain the next layer in the breadth first search.
		set_t nextSet = alloc_set(0);
		//Now propagate from each member of the workingSet.
		for(i = 0; i<size_set(workingSet); i++){

			int currentGroup = workingSet[i];
			set_t currentNeighbors = (set_t) mapsto_set(neighborMatrix, currentGroup);

			///SPECIAL CASE: it is possible that currentNeigbors = NULL because there are no neighbors.
			///This is a special case that can only happen to the exterior empty region
			///when there are no interior empty regions - it still has to have empty cubes
			///but it doesnt necessarily have adjacent empty regions.  Whereas any other
			///empty region would have to have an adjacent empty region in order to ever
			///be referred to in the neighbormatrix. 
			if( currentNeighbors == NULL ){ break; }
			
			//the currentNeighbors are the one sthat we propagate to from currentGroup.
			for(j = 0; j<size_set(currentNeighbors); j++){

				int adjacentGroup = currentNeighbors[j];
				// only propagate if we have not assigned a flag yet.
				if(flags[ adjacentGroup ] != LATTICE_NOT_SET){ continue; }

				// get the box that stores the index to the connectingCube
				int * box = (int *) mapsto_set(currentNeighbors, adjacentGroup);

				// set the connecting cube to what's in the box.
				int connectingCube = box[0];

				//Count the number of intersections
				int numIntersects = cubetoCubeSegmentTest( currentGroup, adjacentGroup, connectingCube );

				//if the number of intersections is even
				if( numIntersects % 2 == 0 ){ flags[adjacentGroup] = flags[currentGroup]; }
				else{
					if( flags[currentGroup] == LATTICE_EXTERIOR ){ flags[adjacentGroup] = LATTICE_INTERIOR; }
					if( flags[currentGroup] == LATTICE_INTERIOR ){ flags[adjacentGroup] = LATTICE_EXTERIOR; }
				}
				// now that we have set a new group, add it to the next set of working groups.
				nextSet = put_set(nextSet, adjacentGroup);
			}
		}

		// now, free the working set and re-allocate,
		// since we want an empty set in the next loop.
		free_set(workingSet);
		workingSet = nextSet;
	}
	free_set(workingSet);

	//now determine if we must use the contingency for the remaining sets
	bool needContingency = false;
	for(i = 0; i<size_set( emptyGroups ); i++){
		if( flags[i] == LATTICE_NOT_SET ){
			needContingency = true;
			break;
		}
	}
	// run contingency testing if necessary.
	if( needContingency ){ flagRemainingSets(); }

	////////////////////////////////////////
//	printf("STATE [populateFlags()    ]: Done.\n");
	////////////////////////////////////////
}






/////////////////////////////////////////////////////////////////////////////////////
// Cube sign generation step 3 (contingency - not always called)
// This function is called only when fillNeighboringSets() is unable to assign a flag to all
// connected components of empty cubes based on neighborhood adjacency.  This function runs
// segment testing between randomly selected cubes within the unflagged group and the
// cubes of the exterior group, which is exterior by definition.
// Note: this should only run on cases of very unususal design.  Unlikely on bio data.
/////////////////////////////////////////////////////////////////////////////////////
void Lattice::flagRemainingSets()
{
	int i = 0;
	int j = 0;

	// first, find the cubeSet that contains 0, 
	// the cube that is empty by definition (because it is the exterior)
	set_t exteriorGroup;
	int containsZero = -1;
	for(i = 0; i<size_set(emptyGroups); i++){
		set_t thisSet = (set_t) mapsto_set(emptyGroups, emptyGroups[i]);
		if( contains_set(thisSet, 0) ){
			containsZero = i;
			exteriorGroup = thisSet;
			break;
		}
	}

	// set this again here, for debugging purposes.
	// this is true by definition.
	if( flags[containsZero] != LATTICE_EXTERIOR ){
		printf("WARNING [flagRemainingSets()  ]: Exterior Group containing the zero cube is not flagged!\n");
		flags[containsZero] = LATTICE_EXTERIOR;	
	}

	////////////////////////////////////////
	printf("WARNING [flagRemainingSets()  ]: Running Contingency testing to assign in/out state... ");
	////////////////////////////////////////
	
	// go through the flags, looking for ones that are not set.	
	for(i = 0; i<size_set(emptyGroups); i++){
		stdoutProgressBar(i, size_set(emptyGroups));
		if( flags[i] != LATTICE_NOT_SET ){ continue; }

	//	printf("EMPTY GROUP %i is not set\n", i);
		set_t unflaggedGroup = (set_t) mapsto_set(emptyGroups, emptyGroups[i]);

		// Make a copy of the unflagged group
		set_t temp = alloc_set(0);	//this will be a copy of unflaggedGroup.
		for(j = 0; j<size_set(unflaggedGroup); j++){
			temp = put_set(temp, unflaggedGroup[j]);
		}

		// fill the unflagged Cubes by randomly removing elements from the copy
		// of the unflagged group.
		set_t unflaggedCubesToTest = alloc_set(0);
		for(j = 0; (j<size_set(unflaggedGroup) && j<NUMBER_OF_CONTINGENCY_TESTS); j++){
			if( size_set(temp)==0 ){ break; }
			int randomVal = random() % size_set(temp);
			unflaggedCubesToTest = put_set(unflaggedCubesToTest, temp[randomVal]);
			remove_set(temp, temp[randomVal]);
		}

		// remake the temp set, and now make it a copy of exteriorGroup.
		// Note that exteriorGroup is actually a member of emptyGroups and must NOT
		// be modified
		free_set(temp); temp = alloc_set(0);
		for(j = 0; j<size_set(exteriorGroup); j++){
			temp = put_set(temp, exteriorGroup[j]);
		}

		// fill the exterior Cubes by randomly removing elements from the copy
		// of the exterior group.
		set_t exteriorCubesToTest = alloc_set(0);
		for(j = 0; (j<size_set(exteriorGroup) && j<NUMBER_OF_CONTINGENCY_TESTS); j++){
			if( size_set(temp)==0 ){ break; }
			int randomVal = random() % size_set(temp);
			exteriorCubesToTest = put_set(exteriorCubesToTest, temp[randomVal]);
			remove_set(temp, temp[randomVal]);
		}

		// now, establish the number of tests to be completed, in case we did not
		// get as many as NUMBER_OF_CONTINGENCY_TESTS
		int numRuns = size_set(unflaggedCubesToTest);
		if( numRuns > size_set(exteriorCubesToTest) ){ numRuns = size_set(exteriorCubesToTest); }

		// run the segment tests.
		int votesInside = 0;
		int votesOutside = 0;
		for(j = 0; j<numRuns; j++){
			int testCube = unflaggedCubesToTest[j];
			int extCube = exteriorCubesToTest[j];
			int output = cubetoCubeSegmentTest( testCube, extCube );
			//printf("NUMBER OF INTERSECTIONS TO EXTERIOR: %i\n", output);
			if( output % 2 == 0 ){ votesOutside++; }	//rememebr extCube is exterior.
			else{ votesInside++; }
		}

		// here we have to vote if it is interior or exterior.
		double val = ((double)votesInside) / ((double)(votesOutside+votesInside));
		if( val>.5 ){ flags[i] = LATTICE_INTERIOR; }
		else{ flags[i] = LATTICE_EXTERIOR; }

		// clean up
		free_set(temp);
		free_set(unflaggedCubesToTest);
		free_set(exteriorCubesToTest);
	}
	// note that we DO NOT free this because it is a 
	// group from emptyGroups, which is used (and freed) later.
	//free_set(exteriorGroup);
}







/////////////////////////////////////////////////////////////////////////////////////
// This function aggregates all the point values, either from caching, cube-flagging, or from
// segment testing.  Then it calls getAverageIntersectionPt() and computes intersection
// points for every edge along two different point values, on a nonempty cube, and stores it
// so that it never needs to be computed again.  Then, on each nonempty cube, it computes
// triangles.  Then it initializes the output surface, and returns.
/////////////////////////////////////////////////////////////////////////////////////
///Cube model used: 
///
///   3---7            ^
///   |\  |\       \   |
///   | 1---5       +y |
///   | | | |        \ z+
///   2-|-6 |         \|  
///    \|  \|          0--x+-->    
///     0---4
///
SurfaceObject * Lattice::finalizeSurface()
{
	////////////////////////////////////////
	printf("STATE [finalizeSurface()    ]: Completing leftover analysis and generating output geometry... ");
	////////////////////////////////////////

	int i = 0;
	int j = 0;

	CubeTable * table = new CubeTable();

	//the points are stored in edgePoints, so store the triangles here.
	set_t triList = alloc_set(SP_MAP);		//arrays of 3 ints, for triangle topology

	int * edgeHolder = new int[12];
	for(i = 0; i<size_set(nonempty); i++){
		stdoutProgressBar(i, size_set(nonempty));

		// Clear the edge holder, so that there is no erroneous data.
		for(j = 0; j<12; j++){ edgeHolder[j]=-1; }

		// pick up the next cube index, for work in thsi loop.
		int cubeIndex = nonempty[i];

		// first, get the cube's indices.	
		int cjkRemainder = cubeIndex % (ydim*zdim);
		int ci = (int) ( ((double) (cubeIndex)) / ((double) (ydim*zdim)  ) );
		int cj = (int) ( ((double) (cjkRemainder)) / ((double) zdim) );
		int ck = cjkRemainder % (zdim);

		// Name the corner indices							///get the corner status for this cube. 
		int c0 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+0) ;      int p0 = isInside( c0 );
		int c1 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+1) ;      int p1 = isInside( c1 );
		int c2 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+0) ;      int p2 = isInside( c2 );
		int c3 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+1) ;      int p3 = isInside( c3 );
		int c4 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+0) ;      int p4 = isInside( c4 );
		int c5 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+1) ;      int p5 = isInside( c5 );
		int c6 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+0) ;      int p6 = isInside( c6 );
		int c7 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+1) ;      int p7 = isInside( c7 );

		// If the corners have different interior/exterior settings, and an edge point has 
		// not been generated in the past, then store an edge point.  This works because "getEdgePoint()" returns -1 if
		// there is no point generated for that edge yet.  This ensures that if we have generated that edge point on a 
		// previous iteration of the loop, we can find it.
		// In the case when we need to store a new edgePoint, storeEdgepoint handles that. 
		// getAverageIntersectionPt() gets the actual geometry of the point for us, and it is passed directly into storage
		if(p0!=p1){ edgeHolder[ 0] = getEdgePoint(c0,c1); if( edgeHolder[ 0] == -1){edgeHolder[ 0] = storeEdgePoint(c0, c1, getAverageIntersectionPt(c0, c1) ); } }
		if(p0!=p2){ edgeHolder[ 1] = getEdgePoint(c0,c2); if( edgeHolder[ 1] == -1){edgeHolder[ 1] = storeEdgePoint(c0, c2, getAverageIntersectionPt(c0, c2) ); } }
		if(p0!=p4){ edgeHolder[ 4] = getEdgePoint(c0,c4); if( edgeHolder[ 4] == -1){edgeHolder[ 4] = storeEdgePoint(c0, c4, getAverageIntersectionPt(c0, c4) ); } }
		if(p4!=p5){ edgeHolder[ 8] = getEdgePoint(c4,c5); if( edgeHolder[ 8] == -1){edgeHolder[ 8] = storeEdgePoint(c4, c5, getAverageIntersectionPt(c4, c5) ); } }
		if(p4!=p6){ edgeHolder[ 9] = getEdgePoint(c4,c6); if( edgeHolder[ 9] == -1){edgeHolder[ 9] = storeEdgePoint(c4, c6, getAverageIntersectionPt(c4, c6) ); } }
		if(p2!=p3){ edgeHolder[ 3] = getEdgePoint(c2,c3); if( edgeHolder[ 3] == -1){edgeHolder[ 3] = storeEdgePoint(c2, c3, getAverageIntersectionPt(c2, c3) ); } }
		if(p2!=p6){ edgeHolder[ 6] = getEdgePoint(c2,c6); if( edgeHolder[ 6] == -1){edgeHolder[ 6] = storeEdgePoint(c2, c6, getAverageIntersectionPt(c2, c6) ); } }
		if(p1!=p3){ edgeHolder[ 2] = getEdgePoint(c1,c3); if( edgeHolder[ 2] == -1){edgeHolder[ 2] = storeEdgePoint(c1, c3, getAverageIntersectionPt(c1, c3) ); } }
		if(p1!=p5){ edgeHolder[ 5] = getEdgePoint(c1,c5); if( edgeHolder[ 5] == -1){edgeHolder[ 5] = storeEdgePoint(c1, c5, getAverageIntersectionPt(c1, c5) ); } }
		if(p6!=p7){ edgeHolder[11] = getEdgePoint(c6,c7); if( edgeHolder[11] == -1){edgeHolder[11] = storeEdgePoint(c6, c7, getAverageIntersectionPt(c6, c7) ); } }
		if(p3!=p7){ edgeHolder[ 7] = getEdgePoint(c3,c7); if( edgeHolder[ 7] == -1){edgeHolder[ 7] = storeEdgePoint(c3, c7, getAverageIntersectionPt(c3, c7) ); } }
		if(p5!=p7){ edgeHolder[10] = getEdgePoint(c5,c7); if( edgeHolder[10] == -1){edgeHolder[10] = storeEdgePoint(c5, c7, getAverageIntersectionPt(c5, c7) ); } }

		// figure out what case we will use in the cubetable.  The cubetable is assignd
		// exponentially, so map the 0-1 interior/exterior status to separate powers of 2.
		int cubeTableIndex = p0 + 2*p1 + 4*p2 + 8*p3 + 16*p4 + 32*p5 + 64*p6 + 128*p7;

		// The indexing gives us coordinates into a matrix in the table, which is how we use
		// the lookup table.
		int numPolysToAdd = table->sizes[cubeTableIndex];
		int * polyList = table->polys[cubeTableIndex];
		
		// Now, iterate through the number of polys we will add.  Now we add the triangles
		for(j = 0; j<numPolysToAdd; j++){
			// allocate the triangle
			int * newTriangle = new int[3];

			// store the indices of the triangle.  Note that since we are dereferencing into
			// the edgeHolder, this is mapping to specific point indices for points we have stored.			
			// reversed order here to reverse norms - original polyList is reversed.
			newTriangle[0] = edgeHolder[ polyList[3*j+0] ];
			newTriangle[1] = edgeHolder[ polyList[3*j+2] ];
			newTriangle[2] = edgeHolder[ polyList[3*j+1] ];

			// store the triangle topology into the triList.  This is how we will save it.
			triList = associate_set(triList, size_set(triList), newTriangle);
		}

		// nothing was allocated, so nothign to delete.
	}

	// Compute statistics
	int numInside = 0;
	int numOutside = 0;
	for(i = 0; i<size_set(points); i++){
		int * test = (int *) mapsto_set(points, points[i]);
		if(test[0] == LATTICE_INTERIOR){ numInside++; }
		if(test[0] == LATTICE_EXTERIOR){ numOutside++; }
	}

	////////////////////////////////////////
	printf("INFO  [Total Interior Points]: %i \n", numInside);
	printf("INFO  [Total Exterior Points]: %i \n", numOutside);
	printf("INFO  [# Geo. Tested Points ]: %i \n", numGeometricChecks);
	printf("INFO  [Total Pts processed  ]: %i \n", size_set(points));
	printf("INFO  [finalizeSurface()    ]: Generating Surface... \n");
	////////////////////////////////////////

	// now generate the resulting surface.
	SurfaceObject * result = new SurfaceObject(edgePoints, triList);

	////////////////////////////////////////
	printf("INFO  [finalizeSurface()    ]: Done generating surface. \n");
	////////////////////////////////////////

	// clean up the table and edgeHOlder that we made, and the triList
	// these data are used by the surfaceObject constructor, but not deleted
	// so we dont need them anymore.
	delete(table);
	delete[](edgeHolder);
	for(i = 0; i<size_set(triList); i++){
		int * temp = (int *) mapsto_set(triList, triList[i]);
		delete[](temp);
	}
	free_set(triList);

	// return the surfaceObject
	return result;
}


















// ##########################################################################################################################
// ###     ##      ##     ##  #####  #      #      #     ##    ##     ######    ##  ##  #     ##     ###    ##     ##      ##
// ##  ###  #  #####  ###  #   ###   #  #######  ###  ##  ##  ##  ###  ####  ##  #  ##  #  ##  #  ##  #  ##  #  ##  ###  ####
// ##  ###  #  #####  ###  #    #    #  #######  ###  ##  ##  ##  #########  #####  ##  #  ##  #  ##  #  ##  #  ##  ###  ####
// ##  ######    ###  ###  #  #   #  #    #####  ###  ##  ##  ##  ##########  ####  ##  #  ##  #  ##  #  ##  #  ##  ###  ####
// ##  ######  #####  ###  #  ## ##  #  #######  ###    ####  ##  ############  ##  ##  #     ##     ##  ##  #    #####  ####
// ##  ##   #  #####  ###  #  #####  #  #######  ###  #  ###  ##  #############  #  ##  #  #####  #####  ##  #  #  ####  ####
// ##  ###  #  #####  ###  #  #####  #  #######  ###  ##  ##  ##  #############  #  ##  #  #####  #####  ##  #  ##  ###  ####
// ##  ###  #  #####  ###  #  #####  #  #######  ###  ##  ##  ##  ###  ####  ##  #  ##  #  #####  #####  ##  #  ##  ###  ####
// ###     ##      ##     ##  #####  #      ###  ###  ##  #    ##     ######    ###    ##  #####  ######    ##  ##  ###  ####
// ##########################################################################################################################

//Heavy Geometry Support Functions//////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Wrapper function for interior/exterior testing.
// First calls getPointStatus, to check info based on the hash and empty cubes.
// This should get the vast majority of points.
// Then calls testShortRay, to do a segment test if the point is totally internal
// Note that this MUST return an answer, and cannot return anything except 
// 1, for inside OR 0 for outside
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::isInside( int pointIndex )
{
	// We perform this test in case the pointStatus can be known from an empty cube.
	// this involves testing 8 cubes, but its much faster than a geometric test.
	// Note that this function returns the point status from cache first, if its cached.
	int firstTest = getPointStatus(pointIndex);

	// If we get data out from the cubes, return it.	
	if( firstTest == LATTICE_INTERIOR ){ return 1; }
	else if( firstTest == LATTICE_EXTERIOR ){ return 0; }

	// If none of the lattice cubes are in the cached, then we need to perform a geometric
	// test.  That is handled below:
	else{
		// When we test the short ray, we find nearby empty cubes and count the number of
		// surface intersections between this point and a point in the center of the cube.
		// printf("RESORTING TO GEOMETRIC TESTS\n");
		double secondTest = testShortRay(pointIndex);
//		if( secondTest > .8 ){ 
		if( secondTest > .5 ){ 
			points = associate_set(points, pointIndex, point_inside);
			numGeometricChecks++;
			return 1; 
		}
//		if( secondTest < .2 ){
		if( secondTest <= .5 ){
			points = associate_set(points, pointIndex, point_outside);
			numGeometricChecks++;
			return 0; 
		}
		else{
			printf("DEBUG: SHORT RAY TESTING FAILED: %f\n", secondTest);
			
			double * pt = getPoint(pointIndex);
			
			printf("FAILED POINT COORDINATES: %f %f %f\n", pt[0], pt[1], pt[2] );
			
			exit(1);
		}
	}
}





/////////////////////////////////////////////////////////////////////////////////////
// shortRay Test:
// This function finds nearby empty cubes, generates a segment to those cubes, and 
// then runs countIntersectionPts to count the number of intersections between the point
// and the cube.  The result is a poll of the number of local cubes.
/////////////////////////////////////////////////////////////////////////////////////
double Lattice::testShortRay(int pointIndex)
{
	int i = 0;

	double * pt = getPoint(pointIndex);
	
	//get the cubes local to the point.
	set_t localCubes = getCubesAdjacentToPoint(pointIndex, LATTICE_RADIUS_TO_TEST);
	//get the cubes adjacent to the point.
	set_t adjCubes = getCubesAdjacentToPoint(pointIndex, 1);
	
	//go through the local indices
	int votesInside = 0;
	int votesOutside = 0;
	
	int debug_numCubesTested = 0;
	for(i = 0; i<size_set(localCubes); i++){
		//if the cube is not empty, don't test it.
		int cubeIndex = localCubes[i];
		if( contains_set(nonempty, cubeIndex) ){ continue; }
		if( contains_set(adjCubes, cubeIndex) ){ continue; }
	
		debug_numCubesTested++;
	
		//if it is empty, then perform a segment test between the point and the cube's center.
		//first, get the cube's indices.	
		int cjkRemainder = cubeIndex % (ydim*zdim);
		int ci = (int) ( ((double) (cubeIndex)) / ((double) (ydim*zdim)) );
		int cj = (int) ( ((double) (cjkRemainder)) / ((double) zdim) );
		int ck = cjkRemainder % (zdim);

		//now get the cube's center's coordinates, with a tiny nudge to avoid some geometric issues
		//(notably going precisely diagonally)	
		double xnudge = ( ( (double) (random()%1000) ) / 1000.0 ) * (.1*res);
		double ynudge = ( ( (double) (random()%1000) ) / 1000.0 ) * (.1*res);
		double znudge = ( ( (double) (random()%1000) ) / 1000.0 ) * (.1*res);

		double * cubeCenter = new double[3];
		cubeCenter[0] = xneg + (ci+.5)*res + xnudge;
		cubeCenter[1] = yneg + (cj+.5)*res + ynudge;
		cubeCenter[2] = zneg + (ck+.5)*res + znudge;

		//now get the cube's flag, since it is empty.
		int cubeFlag = getEmptyCubeStatus(cubeIndex);

		//now run the segment test.
		int numIntersects = countIntersectionpts(pt, cubeCenter);
		
		//decide interior/exterior based on number of intersects and flag status.
		if( (numIntersects % 2 == 0) ){
			if( cubeFlag == LATTICE_EXTERIOR ){
				votesOutside++;
			}
			else if( cubeFlag == LATTICE_INTERIOR ){ 
				votesInside++;
			}
		}
		else if( (numIntersects % 2 == 1) ){
			if( cubeFlag == LATTICE_EXTERIOR ){
				votesInside++;
			}
			else if( cubeFlag == LATTICE_INTERIOR ){
				votesOutside++;
			}
		}
		
		delete[](cubeCenter);
	}
	delete[](pt);
	free_set(localCubes);
	free_set(adjCubes);

	double result = ((double)votesInside) / ((double)(votesOutside+votesInside));

	return result;
}

		



/////////////////////////////////////////////////////////////////////////////////////
// Segment testing:
// This function takes a segment and tests how many times it nontrivially intersects
// the input surface.  This eliminates double-intersections where it strikes exactly
// at the edge between two triangles with an epsilon test - thus, the points must
// be apart from each other.
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::countIntersectionpts(double * a, double * b)
{
	//debug
	int i = 0;
	int j = 0;

	///find the cubes that the segment intersects.
	//set_t cubeList = alloc_set(0);
	//cubeList = getSegmentCubes(a, b, cubeList);
	set_t cubeList = getSegmentCubes(a, b);

	//accumulate the triangles into a set.
	set_t triList = alloc_set(0);
	//go through each cube
	for(i = 0; i<size_set(cubeList); i++){
		// get the index of the cube
		int myCube = cubeList[i];
		//if it's an empty cube, ignore.
		if( !contains_set(nonempty, myCube) ){ 
			//printf("CUBE [%i] IS EMPTY, REJECTED\n", myCube);
			continue; 
		}
		//get the triangle indices in the cube
		set_t cubeTriangles = (set_t) mapsto_set(nonempty, myCube);
		//add all the indices into the list of triangles that we will analyze.
		for(j = 0; j<size_set(cubeTriangles); j++){
			//printf("TRILIST INSERTED: %i\n", cubeTriangles[j]);
			triList = put_set(triList, cubeTriangles[j]);
		}
	}

	//now we begin to count intersections.
	set_t iPoints = alloc_set(SP_MAP);
	for(i = 0; i<size_set(triList); i++){
		//get the three points of the triangle t0 t1 t2
		double * tri = surf->getTriangle(triList[i]); 
		//separate the triangle into individual points
		double * t0 = new double[3];	t0[0] = tri[0];	t0[1] = tri[1];	t0[2] = tri[2];
		double * t1 = new double[3];	t1[0] = tri[3];	t1[1] = tri[4];	t1[2] = tri[5];
		double * t2 = new double[3];	t2[0] = tri[6];	t2[1] = tri[7];	t2[2] = tri[8];

		//compute the intersection
		double * intPt = intersect_Seg_Triangle( a, b, t0, t1, t2);
		
		//If the intersection point exists (i.e. not NULL)
		if( intPt != NULL ){
			//go through the points already stored
			for(j = 0; j<size_set(iPoints); j++){
				//and check to see if any are basically identical
				double * oldPt = (double *) mapsto_set(iPoints, j);
				double dist = vectorSize(intPt, oldPt);
				//if they are, delete the intersection point, set it to NULL, and stop checking.
				if(dist < .00001){
					delete[](intPt);
					intPt = NULL;
					break;
				}
			}
			//Since we have not set it to NULL, it must be new, so add it.
			if(intPt != NULL){
				iPoints = associate_set(iPoints, size_set(iPoints), intPt);
			}
		}
		//clear off the temporary data.
		delete[](tri);
		delete[](t0);	delete[](t1);	delete[](t2);
	}
	
	//count the number of unique intersections
	int result = size_set(iPoints);

	//free the temp data used in this computation
	free_set(cubeList);
	free_set(triList);
	for(i = 0; i<size_set(iPoints); i++){
		double * oldPt = (double *) mapsto_set(iPoints, i);
		delete[](oldPt);
	}
	free_set(iPoints);
	
	return result;
}




/////////////////////////////////////////////////////////////////////////////////////
// Segment testing:
// This function counts the sections along a segment between the centers of two
// cubes (e.g. cubetoCubeSegmentTest( testCube, extCube ) ).  A wrapper for the
// function offers the functionality where members of two groups are specified
// as neighbors to a specific cube, and all ordered pairs of cubes, one from one
// group (e.g. currentGroup) and one from the other (e.g. adjacentGroup) are tested
// using cubetoCubeSegmentTest( testCube, extCube ).
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::cubetoCubeSegmentTest( int testCube, int extCube )
{
//	printf("BBQ\n");
	///find the indices of testCube, ijk.
	int jkRemainder = testCube % (ydim*zdim);
	int i = (int) ( ((double) (testCube)) / ((double) (ydim*zdim))  );
	int j = (int) ( ((double) (jkRemainder)) / ((double) (zdim)) );
	int k = jkRemainder % zdim;
	
	//find coordinates for the center of the cube.
	double * a = new double[3];
	a[0] = xneg + (i+.5)*(res);	a[1] = yneg + (j+.5)*(res);	a[2] = zneg + (k+.5)*(res);
	
	double * nudge = new double[3];
	nudge[0] = (random()%1000)/ (double) 1000;
	nudge[1] = (random()%1000)/ (double) 1000;
	nudge[2] = (random()%1000)/ (double) 1000;
	double * norm = normalizeVector(nudge);
	norm[0] *= .1*res;	norm[1] *= .1*res;	norm[2] *= .1*res;	
	
	//nudge the a vector so it's not perfectly centered.
	a[0] += norm[0];	a[1] += norm[1];	a[2] += norm[2];

	///find the indices of extCube, ijk.
	jkRemainder = extCube % (ydim*zdim);
	i = (int) ( ((double) (extCube)) / ((double) (ydim*zdim))  );
	j = (int) ( ((double) (jkRemainder)) / ((double) (zdim)) );
	k = jkRemainder % zdim;
	
	//find coordinates for the center of the cube.
	double * b = new double[3];
	b[0] = xneg + (i+.5)*(res);	b[1] = yneg + (j+.5)*(res);	b[2] = zneg + (k+.5)*(res);

	nudge[0] = (random()%1000)/ (double) 1000;
	nudge[1] = (random()%1000)/ (double) 1000;
	nudge[2] = (random()%1000)/ (double) 1000;
	delete[](norm);norm = normalizeVector(nudge);
	norm[0] *= .1*res;	norm[1] *= .1*res;	norm[2] *= .1*res;	
	
	//nudge the b vector so it's not perfectly centered.
	b[0] += norm[0];	b[1] += norm[1];	b[2] += norm[2];

	int result = countIntersectionpts(a, b);

	delete[](a);
	delete[](b);
	delete[](nudge);
	delete[](norm);

	return result;
}





/////////////////////////////////////////////////////////////////////////////////////
// Segment testing:
// This function counts the sections along a segment between the centers of two
// cubes (e.g. cubetoCubeSegmentTest( testCube, extCube ) ).  A wrapper for the
// function offers the functionality where members of two groups are specified
// as neighbors to a specific cube, and all ordered pairs of cubes, one from one
// group (e.g. currentGroup) and one from the other (e.g. adjacentGroup) are tested
// using cubetoCubeSegmentTest( testCube, extCube ).
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::cubetoCubeSegmentTest( int currentGroup, int adjacentGroup, int connectingCube )
{
	int i = 0;
	int j = 0;

	//first, get the neighbors to the connectingCube.
	set_t neighbors = getNeighborCubes( connectingCube, nonempty, NULL );

	//we will store the neighbors, based on their groups, here
	set_t group1 = alloc_set(0); //currentGroup
	set_t group2 = alloc_set(0); //adjacentGroup

	//now classify neighboring cubes into one of these two groups, or into neither.
	set_t firstGroup = (set_t) mapsto_set(emptyGroups, currentGroup);
	set_t secndGroup = (set_t) mapsto_set(emptyGroups, adjacentGroup);
	for(i = 0; i<size_set(neighbors); i++){
		int cubeIndex = neighbors[i];
		//these SHOULD be mutually exclusive.  It would be super bad if they were not.
		if( contains_set(firstGroup, cubeIndex) ){ group1 = put_set(group1, cubeIndex); }
		if( contains_set(secndGroup, cubeIndex) ){ group2 = put_set(group2, cubeIndex); }
	}

	int numTests = 0;
	int votesEven = 0;
	int votesOdd = 0;
	bool done = false;
	for(i = 0; i<size_set(group1); i++){
		for(j = 0; j<size_set(group2); j++){
			if( numTests > NUMBER_OF_CROSS_SET_TESTS ){ done = true; }

			int tempVal = cubetoCubeSegmentTest( group1[i], group2[j] );
			if( tempVal % 2 == 0 ){ votesEven++; }
			else{ votesOdd++; }

			numTests++;
		}
		if(done){break;}
	}

	free_set(group1);
	free_set(group2);
	free_set(neighbors);

	int result = -1;
	if( votesEven > votesOdd ){ result = 0; }
	else{ result = 1; }

	return result;
}





/////////////////////////////////////////////////////////////////////////////////////
// Segment testing:
// This function gets the intersections along a single edge in the lattice.
// the input are the indices of the two points along the cube to be tested
// Note that even thouse the points are adjacent to up to four cubes, anything
// that intersects the edge is in all four cubes.  So we get all four, and
// eliminate anything that is not in all four.  Then we test segment-tri intersections
// and get the position of the average intersection.
/////////////////////////////////////////////////////////////////////////////////////
double * Lattice::getAverageIntersectionPt(int pt1, int pt2)
{
	int i = 0;
	int j = 0;

	//get the cubes adjacent to each point.
	set_t cubes1 = getCubesAdjacentToPoint(pt1, 1);
	set_t cubes2 = getCubesAdjacentToPoint(pt2, 1);

	//swap cubes1 for the smaller set, for speed.
	if(size_set(cubes1)>size_set(cubes2)){
		set_t temp = cubes1;
		cubes1 = cubes2;
		cubes2 = temp;
	}

	///this is the list of cubes that touch both pt1 and pt2.
	set_t conservedCubes = alloc_set(0);

	//enumerate the cube neighborhood of the first
	for(i = 0; i<size_set(cubes1); i++){
		int ind = cubes1[i];
		//check to see if the cube is in the second cube neighborhood
		if( contains_set(cubes2, ind) ){
			//if both cube neighborhoods contain the same cube, store it.
			conservedCubes = put_set(conservedCubes, ind);
		}
	}

	//error check - if theres nothing conserved, then they were not adjacent.
	if( size_set(conservedCubes) == 0 ){
		printf("ERROR [getAvgIntersectPt]: corners supplied are not adjacent! exit!\n");
		printf("cubes1:");
		for(i = 0; i<size_set(cubes1); i++){ printf(" %i ", cubes1[i]); }
		printf("\n");
		fflush(stdout);
		
		///find the indices of pt1, abc
		int bcRemainder = pt1 % ((ydim+1)*(zdim+1));
		int a = (int) ( ((double) (pt1)) / ((double) ((ydim+1)*(zdim+1)))  );
		int b = (int) ( ((double) (bcRemainder)) / ((double) (zdim+1)) );
		int c = bcRemainder % (zdim+1);
		printf("cube1 coords: %i %i %i\n", a, b, c);		
		
		printf("cubes2:");
		for(i = 0; i<size_set(cubes2); i++){ printf(" %i ", cubes2[i]); }
		printf("\n");

		///find the indices of pt2, abc
		bcRemainder = pt2 % ((ydim+1)*(zdim+1));
		a = (int) ( ((double) (pt2)) / ((double) ((ydim+1)*(zdim+1)))  );
		b = (int) ( ((double) (bcRemainder)) / ((double) (zdim+1)) );
		c = bcRemainder % (zdim+1);
		printf("cube2 coords: %i %i %i\n", a, b, c);		
		fflush(stdout);

		exit(1);
	}

	//now find the cube with the least triangles, while storing pointers to the triangleSets.
	int numTris = INT_MAX;
	int smallest = -1;
	set_t triangleSets = alloc_set(SP_MAP);
	for(i = 0; i<size_set(conservedCubes); i++){
		set_t triSet = (set_t) mapsto_set(nonempty, conservedCubes[i]);
		triangleSets = associate_set(triangleSets, size_set(triangleSets), triSet);
		//find the cube with least triangles
		if( size_set(triSet) < numTris ){
			numTris = size_set(triSet);
			smallest = i;
		}
	}

	set_t smallestSet = (set_t) mapsto_set(triangleSets, smallest);	
//	//do not free this set, it is used in nonempty.
//	remove_set(triangleSets, smallest);

	//now we will find the set of triangles that are in all the cubes, and store them here.
	set_t conservedTriangles = alloc_set(0);
	for(i = 0; i<size_set(smallestSet); i++){
		//get a triangle
		int thisTri = smallestSet[i];
		bool insert = true;
		//look through the other sets
		for(j = 0; j<size_set(triangleSets); j++){
			if(j==smallest){ continue; }
			set_t triSet = (set_t) mapsto_set(triangleSets, triangleSets[j]);
			//if any set does not contain it, do not insert.
			if( !contains_set(triSet, thisTri) ){ insert = false; break; }
		}
		//if all sets contain it, then insert is still true, so insert it.
		if(insert){
			conservedTriangles = put_set(conservedTriangles, thisTri);
		}
	}

	//Now we have the set of triangles that are in ALL the adjacent cubes. 
	//They must intersect the segment.  so now get all the intersection points.
	double * p1 = getPoint(pt1);
	double * p2 = getPoint(pt2);
	set_t intersectionPts = alloc_set(SP_MAP);
	for(i = 0; i<size_set(conservedTriangles); i++){
		int triIndex = conservedTriangles[i];
		double * tri = surf->getTriangle(triIndex); 
		//separate the triangle into individual points
		double * t0 = new double[3];	t0[0] = tri[0];	t0[1] = tri[1];	t0[2] = tri[2];
		double * t1 = new double[3];	t1[0] = tri[3];	t1[1] = tri[4];	t1[2] = tri[5];
		double * t2 = new double[3];	t2[0] = tri[6];	t2[1] = tri[7];	t2[2] = tri[8];
	
		//get the intersectionPt
		double * thisPoint = intersect_Seg_Triangle( p1, p2, t0, t1, t2);

		//I dont think a NULL should happen, but catch it.
		if( thisPoint != NULL){
			for(j = 0; j<size_set(intersectionPts); j++){
				//and check to see if any are basically identical
				double * oldPt = (double *) mapsto_set(intersectionPts, j);
				double dist = vectorSize(thisPoint, oldPt);
				//if they are, delete the intersection point, set it to NULL, stop checking.
				if(dist < .00001){
					delete[](thisPoint);
					thisPoint = NULL;
					break;
				}
			}
			//Since we have not set it to NULL, it must be new, so add it.
			if(thisPoint != NULL){
				intersectionPts = associate_set(intersectionPts, size_set(intersectionPts), thisPoint);
			}
		}		
		delete[](tri);	delete[](t0);	delete[](t1);	delete[](t2);
	}

	if( size_set(intersectionPts)==0 ){
		double * temp = new double[3];
		temp[0] = (p1[0]+p2[0])/2;
		temp[1] = (p1[1]+p2[1])/2;
		temp[2] = (p1[2]+p2[2])/2;
		intersectionPts = associate_set(intersectionPts, size_set(intersectionPts), temp);
		
	//	printf("ERROR [getAverageIntersectionPts]: no intersection found between these points! exit!\n");
	//	exit(1);
	}
	
	///Now that we have the intersection points, average them.
	double * result = new double[3];
	result[0] = 0;	result[1] = 0;	result[2] = 0;

	for(i = 0; i<size_set(intersectionPts); i++){
		double * myPt = (double *) mapsto_set(intersectionPts, i);
		result[0] += myPt[0];	result[1] += myPt[1];	result[2] += myPt[2];
	}
	result[0] /=size_set(intersectionPts);
	result[1] /=size_set(intersectionPts);
	result[2] /=size_set(intersectionPts);

	//clean up
	free_set(cubes1);
	free_set(cubes2);
	free_set(conservedCubes);
	free_set(triangleSets);
	free_set(conservedTriangles);
	delete[](p1);
	delete[](p2);
	for(i = 0; i<size_set(intersectionPts); i++){
		double * pt = (double *) mapsto_set(intersectionPts, i);
		delete[](pt);
	}
	free_set(intersectionPts);

	return result;	
}


















// ##########################################################################################################################
// ##  ###  #      #  #####     ##      #     #####      #  ##  #  ###  ##    ##      #    ##    ##  ###  ##    #############
// ##  ###  #  #####  #####  ##  #  #####  ##  ####  #####  ##  #   ##  #  ##  ###  ####  ##  ##  #   ##  #  ##  ############
// ##  ###  #  #####  #####  ##  #  #####  ##  ####  #####  ##  #    #  #  #######  ####  ##  ##  #    #  #  ################
// ##       #    ###  #####  ##  #    ###  ##  ####     ##  ##  #  #    #  #######  ####  ##  ##  #  #    ##  ###############
// ##  ###  #  #####  #####     ##  #####    ######  #####  ##  #  ##   #  #######  ####  ##  ##  #  ##   ####  #############
// ##  ###  #  #####  #####  #####  #####  #  #####  #####  ##  #  ###  #  #######  ####  ##  ##  #  ###  #####  ############
// ##  ###  #  #####  #####  #####  #####  ##  ####  #####  ##  #  ###  #  #######  ####  ##  ##  #  ###  #####  ############
// ##  ###  #  #####  #####  #####  #####  ##  ####  #####  ##  #  ###  #  ##  ###  ####  ##  ##  #  ###  #  ##  ############
// ##  ###  #      #      #  #####      #  ##  ####  ######    ##  ###  ##    ####  ###    ##    ##  ###  ##    #############
// ##########################################################################################################################

//Helper Functions/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Helper Function - Given the index of a cube, returns the indices of
// cubes adjacent to this one, or returns NULL if it is out of range.
// The cubes are virtual.  They are conceptually indexed in a grid, with indices
//    i*(ydim*zdim) + j*(zdim) + k;
// - will not include an adjacent cube in set_t "excluded", if this is non-null;
// if outputSet is NULL, allocate a new one.
// Here adjacency applies ONLY TO CUBES THAT SHARE A FACE.
/////////////////////////////////////////////////////////////////////////////////////
set_t Lattice::getAdjacentCubes( int ind, set_t excluded, set_t outputSet )
{
	///if outputSet is NULL, allocate a new one.
	if(outputSet == NULL){ outputSet = alloc_set(0); }
	
	///find the indices of this cube, ijk.
	int jkRemainder = ind % (ydim*zdim);
	int i = (int) ( ((double) (ind)) / ((double) (ydim*zdim))  );
	int j = (int) ( ((double) (jkRemainder)) / ((double) (zdim)) );
	int k = jkRemainder % zdim;

	///check bounds and return NULL if broken.
	if( ind<0 || ind>=(xdim*ydim*zdim) || i<0 || i>=xdim || j<0 || j>=ydim || k<0 || k>=zdim ){
		printf("ERROR [getAdjacentCubes()]: supplied cube index, %i, is out of bounds\n", ind);
		return NULL;
	}

	if( excluded == NULL ){
		if(i != 0)     { outputSet = put_set(outputSet, (i-1)*(ydim*zdim) + (j*zdim) + k); }
		if(i != xdim-1){ outputSet = put_set(outputSet, (i+1)*(ydim*zdim) + (j*zdim) + k); }
		if(j != 0)     { outputSet = put_set(outputSet, i*(ydim*zdim) + ((j-1)*zdim) + k); }
		if(j != ydim-1){ outputSet = put_set(outputSet, i*(ydim*zdim) + ((j+1)*zdim) + k); }
		if(k != 0)     { outputSet = put_set(outputSet, i*(ydim*zdim) + (j*zdim) + (k-1)); }
		if(k != zdim-1){ outputSet = put_set(outputSet, i*(ydim*zdim) + (j*zdim) + (k+1)); }
	}
	else{
		if( i != 0 && !contains_set(excluded, (i-1)*(ydim*zdim) + (j*zdim) + k) ){ 
			outputSet = put_set(outputSet, (i-1)*(ydim*zdim) + (j*zdim) + k); 
		}
		if( i != xdim-1 && !contains_set(excluded, (i+1)*(ydim*zdim) + (j*zdim) + k) ){ 
			outputSet = put_set(outputSet, (i+1)*(ydim*zdim) + (j*zdim) + k);
		}
		if( j != 0 && !contains_set(excluded, i*(ydim*zdim) + ((j-1)*zdim) + k) ){ 
			outputSet = put_set(outputSet, i*(ydim*zdim) + ((j-1)*zdim) + k);
		}
		if( j != ydim-1 && !contains_set(excluded, i*(ydim*zdim) + ((j+1)*zdim) + k) ){
			outputSet = put_set(outputSet, i*(ydim*zdim) + ((j+1)*zdim) + k);
		}
		if( k != 0 && !contains_set(excluded, i*(ydim*zdim) + (j*zdim) + (k-1)) ){
			outputSet = put_set(outputSet, i*(ydim*zdim) + (j*zdim) + (k-1));
		}
		if( k != zdim-1 && !contains_set(excluded, i*(ydim*zdim) + (j*zdim) + (k+1)) ){
			outputSet = put_set(outputSet, i*(ydim*zdim) + (j*zdim) + (k+1));
		}
	}

	return outputSet;

	//This code is saved in case a version of this function is needed that outputs an array.
/*	int * result = new int[6];
	if(i == 0){ result[0] = -1; }
	else{ result[0] = (i-1)*(ydim*zdim) + (j*zdim) + k; }
	if(i == xdim-1){ result[1] = -1; }
	else{ result[1] = (i+1)*(ydim*zdim) + (j*zdim) + k; }
	if(j == 0){ result[2] = -1; }
	else{ result[2] = i*(ydim*zdim) + ((j-1)*zdim) + k; }
	if(j == ydim-1){ result[3] = -1; }
	else{ result[3] = i*(ydim*zdim) + ((j+1)*zdim) + k; }
	if(k == 0){ result[4] = -1; }
	else{ result[4] = i*(ydim*zdim) + (j*zdim) + (k-1); }
	if(k == zdim-1){ result[5] = -1; }
	else{ result[5] = i*(ydim*zdim) + (j*zdim) + (k+1); }

	if (excluded != NULL){
		if( result[0] != -1 && contains_set(excluded, result[0]) ){ result[0] = -1; }
		if( result[1] != -1 && contains_set(excluded, result[1]) ){ result[1] = -1; }
		if( result[2] != -1 && contains_set(excluded, result[2]) ){ result[2] = -1; }
		if( result[3] != -1 && contains_set(excluded, result[3]) ){ result[3] = -1; }
		if( result[4] != -1 && contains_set(excluded, result[4]) ){ result[4] = -1; }
		if( result[5] != -1 && contains_set(excluded, result[5]) ){ result[5] = -1; }
	}
*/

}






/////////////////////////////////////////////////////////////////////////////////////
// Helper function - this gives the set of cubes that shares a face or a
// point with the cube provided (by it's index), except for any excluded
// cube.  Results are put into outputSet, unless outputSet is NULL, in
// which case a new set is allocated and returned.
/////////////////////////////////////////////////////////////////////////////////////
set_t Lattice::getNeighborCubes( int index, set_t excluded, set_t outputSet )
{
	///if outputSet is NULL, allocate a new one.
	if(outputSet == NULL){ outputSet = alloc_set(0); }
	
	///find the indices of this cube, ijk.
	int jkRemainder = index % (ydim*zdim);
	int i = (int) ( ((double) (index)) / ((double) (ydim*zdim))  );
	int j = (int) ( ((double) (jkRemainder)) / ((double) (zdim)) );
	int k = jkRemainder % zdim;
	
	///find the nearby range.
	int ilo = i-1;	int ihi = i+1;		if(ilo<0){ ilo = i; }	if(ihi>=xdim){ ihi=i; }
	int jlo = j-1;	int jhi = j+1; 	if(jlo<0){ jlo = i; }	if(ihi>=xdim){ jhi=i; }
	int klo = k-1;	int khi = k+1; 	if(klo<0){ klo = i; }	if(ihi>=xdim){ khi=i; }
	
	int a, b, c;
	for(a = ilo; a<=ihi; a++){
		for(b = jlo; b<=jhi; b++){
			for(c = klo; c<=khi; c++){
				int testValue = a*(ydim*zdim) + b*(zdim) + c;
				if( !contains_set(excluded, testValue) && testValue != index){
					outputSet = put_set(outputSet, testValue);
				}
			}
		}
	}
	
	return outputSet;
}






/////////////////////////////////////////////////////////////////////////////////////
// Helper Function - Union find on cube adjacency.
// This performs a rapid union find on the cubes that are indexed in set_t cubes.
// The assumed ajacency rules are applied, except in the case of being adjacent
// to something in excluded - those cannot be part of any union.
// Returns a set of sets, which define the equivalence classes.
/////////////////////////////////////////////////////////////////////////////////////
set_t Lattice::unionFind( set_t excluded )
{
	int i = 0;
	
	set_t result = alloc_set(SP_MAP);
	set_t unused = alloc_set(0);
	set_t wavefront = alloc_set(0);
	set_t tempSet = alloc_set(0);

	///fill the array with unused cubes that are not excluded.
	for(i = 0; i<xdim*ydim*zdim; i++){
		if( !contains_set(excluded, i) ){
			unused = put_set(unused, i);
		}
	}
	int numUnused = size_set(unused);
	
	///Run Union Find until there is nothing left.
	while( size_set(unused) > 0 ){
		set_t newGroup = alloc_set(0);

		wavefront = put_set(wavefront, unused[0]);
		while(size_set(wavefront) > 0){
			///put everything in the wavefront into the newGroup
			for(i = 0; i<size_set(wavefront); i++){
				stdoutProgressBar(numUnused-size_set(unused), numUnused);
				newGroup = put_set(newGroup, wavefront[i]);
				remove_set(unused, wavefront[i]);
			}

			///get all cubes adjacent to the wavefront, except the excluded, add to the tempSet
			for(i = 0; i<size_set(wavefront); i++){
				tempSet = getAdjacentCubes( wavefront[i], excluded, tempSet);
			}

			///free the wavefront for the next wave.			
			free_set(wavefront);
			wavefront = alloc_set(0);

			///copy the tempset into the wavefront, to make the next wave.
			for(i = 0; i<size_set(tempSet); i++){
				//dont put in stuff that is already in the newgroup.
				if( !contains_set(newGroup, tempSet[i]) ){
					wavefront = put_set(wavefront, tempSet[i]);
				}
			}

			///free the tempset for next usage
			free_set(tempSet);
			tempSet = alloc_set(0);
		}

		result = associate_set(result, size_set(result), newGroup);
	}

/*
//	-----------------------------DEBUGGING--------------------------
	///print out the unions.
	printf("TOTAL NUMBER OF GROUPS FOUND: %i\n", size_set(result) );
	int j = 0;
	set_t firstSet = (set_t) mapsto_set(result, 0);
	for(i = 0; i<size_set(result); i++){
		set_t printSet = (set_t) mapsto_set(result, i);
		printf("GROUP %i: %i elements.\n", i, size_set(printSet) );
//		printf("{ ", i);
		if(i>0){
			printf("ASD1\n");
			for(j = 0; j<size_set(printSet); j++){
//				printf("ASD2\n");
				if( contains_set(firstSet, printSet[j]) ){
					printf("WTF SET # %i CONTAINS SOMETHING IN FIRST-SET (%i)\n", j, printSet[j]);
				}
			}
			printf("ASD3\n");
			//printf("%i, ", printSet[j]);
		}
//		printf("}\n");
		printf("ASD4\n");
	}
	printf("ASD\n");
//	-----------------------------DEBUGGING--------------------------
*/

	free_set(tempSet);
	free_set(wavefront);
	free_set(unused);
	
	return result;
}





/////////////////////////////////////////////////////////////////////////////////////
// Helper Function - gets the cubes that this segment intersects.  exact.
// note that the segment will almost always start with point a on a cube corner.
//  THis is called recursively until it gets to the end.
//  IMPORTANT: only works if both endpoints are inside the lattice.
//  Wrong or bugged otherwise!!!
//  ALSO it is assumed that v1 is generally placed, but v2 is DEFINITELY NOT A CORNER.
/////////////////////////////////////////////////////////////////////////////////////
set_t Lattice::getSegmentCubes(double * v1, double * v2)
{
//	int i = 0;
//	int j = 0;
	
	set_t result = alloc_set(0);
	int xind, yind, zind;

//	printf("RUNNING GETSEGMENTCUBES--------------------------------------------------------\n");
//	printf("xneg: %f  yneg: %f  zneg: %f   res: %f\n", xneg, yneg, zneg, res);
	
	double fmodx = fabs(fmod( v1[0]-xneg, res ));
	double fmody = fabs(fmod( v1[1]-yneg, res ));
	double fmodz = fabs(fmod( v1[2]-zneg, res ));

	if( fmodx<.00001 || fmody<.00001 || fmodz<.00001 ){
		//the point is a corner, so nudge closer to v1 so that it's inside the cube that the
		//segment will pass through, and find out what segment it is.
		double * norm = normalizeVector( v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2] );
		double nx = v1[0] + ((res*.01)*norm[0]);
		double ny = v1[1] + ((res*.01)*norm[1]);
		double nz = v1[2] + ((res*.01)*norm[2]);
//		printf("I'm a corner\n");
//		printf("nx: %f  ny: %f  nz: %f\n", nx, ny, nz);

		xind = (int) ((nx-xneg)/res);	//if(xind == xdim){ xind = xdim-1; }//these not necessary
		yind = (int) ((ny-yneg)/res);	//if(yind == ydim){ yind = ydim-1; }//since we are inside
		zind = (int) ((nz-zneg)/res);	//if(zind == zdim){ zind = zdim-1; }//the cube via nudge (nx,ny,nz)
		delete[](norm);
	}
	else{
//		printf("I'm not a corner\n");
		xind = (int) ((v1[0]-xneg)/res);	//if(xind == xdim){ xind = xdim-1; }//these not necessary
		yind = (int) ((v1[1]-yneg)/res);	//if(yind == ydim){ yind = ydim-1; }//because we dont
		zind = (int) ((v1[2]-zneg)/res);	//if(zind == zdim){ zind = zdim-1; }//need a nudge anyway
	}


	///put the first point in.
	result = put_set(result, xind*(ydim*zdim) + yind*zdim + zind );

//	int counter = 0;

//	printf("############################################# STARTING LOOP\n");
//	printf("v1: %f %f %f\n", v1[0], v1[1], v1[2] );
//	printf("v2: %f %f %f\n", v2[0], v2[1], v2[2] );
//	printf("xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f \n", 
//			xneg, xpos, yneg, ypos, zneg, zpos );
//	printf("RESULT: [");
//	for(i = 0; i<size_set(result); i++){
//		printf("%i ", result[i]);
//	}
//	printf("]\n");

	while( 1 ){
//		printf("############################################# resultSize: %i\n", size_set(result) );

//		counter++;
//		printf("counter: %i\n", counter);
		int xbot=0, xtop=0, ybot=0, ytop=0, zbot=0, ztop=0;

		// declare the cube constraints
		double xlo = xneg + ((xind+0)*res);	double xhi = xneg + ((xind+1)*res);
		double ylo = yneg + ((yind+0)*res);	double yhi = yneg + ((yind+1)*res);
		double zlo = zneg + ((zind+0)*res);	double zhi = zneg + ((zind+1)*res);

		// find out what walls it hits, unless we are at a border cube.
		if( xind!=0      ){ xbot = segAAquadCollision(v1, v2, xlo, xlo, ylo, yhi, zlo, zhi, NULL); }
		if( xind!=xdim-1 ){ xtop = segAAquadCollision(v1, v2, xhi, xhi, ylo, yhi, zlo, zhi, NULL); }
		if( yind!=0      ){ ybot = segAAquadCollision(v1, v2, xlo, xhi, ylo, ylo, zlo, zhi, NULL); }
		if( yind!=ydim-1 ){ ytop = segAAquadCollision(v1, v2, xlo, xhi, yhi, yhi, zlo, zhi, NULL); }
		if( zind!=0      ){ zbot = segAAquadCollision(v1, v2, xlo, xhi, ylo, yhi, zlo, zlo, NULL); }
		if( zind!=zdim-1 ){ ztop = segAAquadCollision(v1, v2, xlo, xhi, ylo, yhi, zhi, zhi, NULL); }

//		//debug
//		printf("1) %i Cube sides hit:", xbot+xtop+ybot+ytop+zbot+ztop);
//		if(xbot){ printf(" xbot"); }
//		if(xtop){ printf(" xtop"); }
//		if(ybot){ printf(" ybot"); }
//		if(ytop){ printf(" ytop"); }
//		if(zbot){ printf(" zbot"); }
//		if(ztop){ printf(" ztop"); }
//		printf("\n");
//		printf("xlo %f xhi %f ylo %f yhi %f zlo %f zhi %f \n", xlo, xhi, ylo, yhi, zlo, zhi );
//		printf("xind: %i (xdim: %i) yind: %i (ydim: %i) zind: %i (zdim: %i) \n", xind, xdim, yind, ydim, zind, zdim);


		// eliminate cases where v1 is on a corner.
		if( xbot==1 && fabs(v1[0]-xlo)<.00001 ){ xbot = 0; }
		if( xtop==1 && fabs(v1[0]-xhi)<.00001 ){ xtop = 0; }
		if( ybot==1 && fabs(v1[1]-ylo)<.00001 ){ ybot = 0; }
		if( ytop==1 && fabs(v1[1]-yhi)<.00001 ){ ytop = 0; }
		if( zbot==1 && fabs(v1[2]-zlo)<.00001 ){ zbot = 0; }
		if( ztop==1 && fabs(v1[2]-zhi)<.00001 ){ ztop = 0; }

		//the segment hits the cube in two places, one of which comes from the cube we were just in.
		//that cube is already in the result list - so check for it, and then remove it.
		if(xbot==1 && contains_set(result, (xind-1)*(ydim*zdim) + (yind+0)*zdim + (zind+0) ) ){ xbot = 0; }
		if(xtop==1 && contains_set(result, (xind+1)*(ydim*zdim) + (yind+0)*zdim + (zind+0) ) ){ xtop = 0; }
		if(ybot==1 && contains_set(result, (xind+0)*(ydim*zdim) + (yind-1)*zdim + (zind+0) ) ){ ybot = 0; }
		if(ytop==1 && contains_set(result, (xind+0)*(ydim*zdim) + (yind+1)*zdim + (zind+0) ) ){ ytop = 0; }
		if(zbot==1 && contains_set(result, (xind+0)*(ydim*zdim) + (yind+0)*zdim + (zind-1) ) ){ zbot = 0; }
		if(ztop==1 && contains_set(result, (xind+0)*(ydim*zdim) + (yind+0)*zdim + (zind+1) ) ){ ztop = 0; }

//		//debug
//		printf("2) %i Cube sides hit", xbot+xtop+ybot+ytop+zbot+ztop);
//		if(xbot){ printf(" xbot"); }
//		if(xtop){ printf(" xtop"); }
//		if(ybot){ printf(" ybot"); }
//		if(ytop){ printf(" ytop"); }
//		if(zbot){ printf(" zbot"); }
//		if(ztop){ printf(" ztop"); }
//		printf("\n");

		//error check
		if( xbot+xtop+ybot+ytop+zbot+ztop>2 ){ 
			printf("ERROR [getSegmentCubes]: segment intersected %i cubes! WTF\n", xbot+xtop+ybot+ytop+zbot+ztop);
			if(xbot){ printf("xbot\n"); }
			if(xtop){ printf("xtop\n"); }
			if(ybot){ printf("ybot\n"); }
			if(ytop){ printf("ytop\n"); }
			if(zbot){ printf("zbot\n"); }
			if(ztop){ printf("ztop\n"); }
			
			printf("v1: %f %f %f\n", v1[0], v1[1], v1[2] );
			printf("v2: %f %f %f\n", v2[0], v2[1], v2[2] );
			printf("xlo %f xhi %f ylo %f yhi %f zlo %f zhi %f \n", xlo, xhi, ylo, yhi, zlo, zhi );
			//printf("counter = %i\n");
			
			//if(v1Corner){ printf("V1 CORNER true\n"); }
		}


		//if one of these is 1, then we can go in that dir bc it's not on the wall already.
		//becayse we checked that above)
		int resultSize = size_set(result);
		if( xbot==1 ){ xind = xind-1; result = put_set( result, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( xtop==1 ){ xind = xind+1; result = put_set( result, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( ybot==1 ){ yind = yind-1; result = put_set( result, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( ytop==1 ){ yind = yind+1; result = put_set( result, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( zbot==1 ){ zind = zind-1; result = put_set( result, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( ztop==1 ){ zind = zind+1; result = put_set( result, xind*(ydim*zdim) + yind*zdim + zind ); }

//		printf("CUBE RESULT SET SIZE: %i\n", size_set(result));
//		printf("xlo %f xhi %f ylo %f yhi %f zlo %f zhi %f \n", xlo, xhi, ylo, yhi, zlo, zhi );
//		printf("xind: %i (xdim: %i) yind: %i (ydim: %i) zind: %i (zdim: %i) \n", xind, xdim, yind, ydim, zind, zdim);
//
//		printf("RESULT: [");
//		for(i = 0; i<size_set(result); i++){
//			if( contains_set(nonempty, result[i]) ){ printf("<%i> ", result[i]); }
//			else for(j = 0; j<size_set(emptyGroups); j++){
//				set_t temp = (set_t) mapsto_set(emptyGroups, emptyGroups[j]);
//				if( contains_set(temp, result[i]) ){
//					printf("%i(%i) ", result[i], emptyGroups[j]);
//				}
//			}	
//		}
//		printf("]\n");
//

		if( size_set(result) == resultSize ){ 
			break;
		}
		
	}
	
//	printf("DONE RUNNING GETSEGMENTCUBES---------------------------------------------------\n");
	return result;
}





/////////////////////////////////////////////////////////////////////////////////////
// Helper Function - Finds out if the point in question is inside or outside based
// on assigned cubes.  Beginning with the index of the point, which are indexed as follows,
//    i*((ydim+1)*(zdim+1)) + j*(zdim+1) + k;
// we then find the index of the cubes it is adjacent to, and check their int/ext status.
// returns LATTICE_INTERIOR if the point is inside, LATTICE_EXTERIOR of the point is outside
// and LATTICE_NONEMPTY if the point is surrounded by nonempties.
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::getPointStatus(int pointIndex)
{
	//first check hash.
	if( contains_set(points, pointIndex) ){
		int * state = (int *) mapsto_set(points, pointIndex);
		return state[0];
	}
	
	//if it's not in the hash already, look at adjacent cubes.
	int i = 0;
	int j = 0;
	int result = LATTICE_NONEMPTY;
	set_t adjacentCubes = getCubesAdjacentToPoint(pointIndex, 1);
	
	for(i = 0; i<size_set(adjacentCubes); i++){
		int ind = adjacentCubes[i];
		//if the adjacent cube is nonempty, it cant be used to assign signs. skip.
		if( contains_set(nonempty, ind) ){ continue; }

		//then it must be empty, so find out the empty cube's sign.		
		for(j = 0; j<size_set(emptyGroups); j++){
			set_t emptyGroup = (set_t) mapsto_set(emptyGroups, j);
			///when we find the empty group with the ind, get the flag info and return.
			if( contains_set(emptyGroup, ind) ){
				result = flags[j];
				break;
			}
		}
		
		//if we've found the flag, then we dont have to look at more cubes.
		if( result != LATTICE_NONEMPTY){ break; }
	}
	
	// if we've found a flag, store it.
	if( result == LATTICE_INTERIOR ){ points = associate_set(points, pointIndex, point_inside); }
	if( result == LATTICE_EXTERIOR ){ points = associate_set(points, pointIndex, point_outside); }
	
	free_set(adjacentCubes);
	
	// if we have not found the flag, then the result will still be LATTICE_NOT_SET, since it 
	// is adjacent to only nonempty cubes.
	return result;	
}






/////////////////////////////////////////////////////////////////////////////////////
// Helper function finds the cubes adjacent to a point, given the point index.
/////////////////////////////////////////////////////////////////////////////////////
set_t Lattice::getCubesAdjacentToPoint(int pointIndex, int radius)
{
	set_t result = alloc_set(0);
	
	if(radius == 0){ return result; }
	
	///find the indices of this point, xyz.
	int jkRemainder = pointIndex % ((ydim+1)*(zdim+1));
	int x = (int) ( ((double) (pointIndex)) / ((double) (ydim+1)*(zdim+1))  );
	int y = (int) ( ((double) (jkRemainder)) / ((double) (zdim+1)) );
	int z = jkRemainder % (zdim+1);
	
	///not that these ranges are defined inclusively
	int xlo = x-radius;
	int ylo = y-radius;
	int zlo = z-radius;
	int xhi = x+(radius-1);
	int yhi = y+(radius-1);
	int zhi = z+(radius-1);

	///note that these ranges are thus also bounded inclusively
	if(xlo<0){ xlo = 0; }
	if(ylo<0){ ylo = 0; }
	if(zlo<0){ zlo = 0; }
	if(xhi>=xdim){ xhi = xdim-1; }
	if(yhi>=ydim){ yhi = ydim-1; }
	if(zhi>=zdim){ zhi = zdim-1; }
	
	//thus we must iterate inclusively.
	int i=0, j=0, k=0;
	for(i = xlo; i<=xhi; i++){
		for(j = ylo; j<=yhi; j++){
			for(k = zlo; k<=zhi; k++){
				result = put_set(result, i*(ydim*zdim) + j*zdim + k);
			}
		}
	}

	return result;
}






/////////////////////////////////////////////////////////////////////////////////////
// Helper function gets the double vector position of a point using it's index.
/////////////////////////////////////////////////////////////////////////////////////
double * Lattice::getPoint(int ind)
{
	// find the indices of this point, ijk.
	int jkRemainder = ind % ((ydim+1)*(zdim+1));
	int i = (int) ( ((double) (ind)) / ((double) (ydim+1)*(zdim+1))  );
	int j = (int) ( ((double) (jkRemainder)) / ((double) (zdim+1)) );
	int k = jkRemainder % (zdim+1);
	
	// Instantiate theresult, and fill it with the coordinates.
	double * result = new double[3];
	result[0] = xneg + (i*res);
	result[1] = yneg + (j*res);
	result[2] = zneg + (k*res);

	// Return the final output.s
	return result;
}



/////////////////////////////////////////////////////////////////////////////////////
// Helper function gets the flag assigned to a specific cube
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::getEmptyCubeStatus(int cubeIndex)
{
	if( cubeIndex < 0 || cubeIndex >= xdim*ydim*zdim ){
		printf("ERROR [getCubeFlag]: cubeIndex provided was out of range! exit!\n");
		exit(1);
	}
	
	//first check nonempty, since we can check that fast.
	if( contains_set(nonempty, cubeIndex) ){ 
		printf("OMG WTF THE CUBE YOU GAVE ME IS NOT EMPTY!\n");
		exit(1);
	}
	
	int i = 0;
	int result = 0;
	//next check each group of empties.
	for(i = 0; i<size_set(emptyGroups); i++){
		set_t thisGroup = (set_t) mapsto_set(emptyGroups, i);
		if( contains_set(thisGroup, cubeIndex) ){
			result = flags[i];
		}
	}

	return result;
}




/////////////////////////////////////////////////////////////////////////////////////
// Helper functions get and store edgePoints into edges
// Returns the index into edgePoints for this edge POint.
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::storeEdgePoint(int p0, int p1, double * pt)
{
	///a must always be lower, so we can save space.
	int a = p0;
	int b = p1;
	if(a > b){
		a = p1;
		b = p0;
	}
	
	set_t aSet;
	if( contains_set(edges, a) ){	aSet = (set_t) mapsto_set(edges, a); }
	else{ aSet = alloc_set(SP_MAP); }
	
	if( contains_set(aSet, b) ){
		printf("ERROR [storeEdgePoint]: already contains an edge here!\n");
		exit(1);
	}
	
	int thisIndex = size_set(edgePoints);
	//now that we can inser the point, insert it.
	edgePoints = associate_set(edgePoints, thisIndex, pt);

	//store the index in this array
	int * indexHolder = new int[1];
	indexHolder[0] = thisIndex;

	//store the array in edges
	aSet = associate_set(aSet, b, indexHolder);
	edges = associate_set(edges, a, aSet);
	
	return thisIndex;
}




/////////////////////////////////////////////////////////////////////////////////////
//returns -1 if point pair is not stored.
//returns the index into edgePOints if the pair is stored correctly.
/////////////////////////////////////////////////////////////////////////////////////
int Lattice::getEdgePoint(int p0, int p1)
{
	///a must always be lower, so we can save space.
	int a = p0;
	int b = p1;
	if(a > b){
		a = p1;
		b = p0;
	}

	set_t aSet;
	int result = -1;
	
	//check to see that it exists, if not, return NULL
	if( contains_set(edges, a) ){	
		aSet = (set_t) mapsto_set(edges, a); 
	}
	else{ return -1; }
	
	//check to see if it exists, if not return NULL
	if( !contains_set(aSet, b) ){ return -1; }
	else{
		//if it does, get the box, get the index, and get the double, and return it.
		int * box = (int *) mapsto_set(aSet, b);
		result = box[0];
	}
	
	return result;
}




/////////////////////////////////////////////////////////////////////////////////////
// helper function returns a set containing the indices of the empty groups that
// have cubes adjacent to the input cube.  The input cube is expected to be a *NON-EMPTY*
// cube, and thus finding more than one empty groups on either side HELPS us to find
// NON-empty regions that connect empty regions. This sees if the regions on the other sides
// are differently indexed (e.g. topologically non-identical.)
/////////////////////////////////////////////////////////////////////////////////////
set_t Lattice::getNeighboringEmptyGroups( int cubeIndex )
{
	int i = 0;
	int j = 0;
	set_t result = alloc_set(0);
	
	//first get the neighbors
	set_t nbors = getNeighborCubes( cubeIndex, nonempty, NULL);
	for(i = 0; i<size_set(nbors); i++){
		int myCube = nbors[i];
		for(j = 0; j<size_set(emptyGroups); j++){
			set_t thisGroup = (set_t) mapsto_set(emptyGroups, j);
			if( contains_set(thisGroup, myCube) ){
				result = put_set(result, j);
				break;
			}
		}
	}
	free_set(nbors);	
	
	return result;
}


// #####################################################################################################################################
// #####################################################################################################################################
// ####################    ###########    ######            ##            ##        #####       #####         ##########################
// ####################    ##########      #####            ##            ##        ###           ###         ##########################
// ####################    #########        ########    ##########    ########    ####    #####    ##    ###############################
// ####################    ########    ##    #######    ##########    ########    ####    ###########    ###############################
// ####################    #######    ####    ######    ##########    ########    ####    ###########       ############################
// ####################    #######    ####    ######    ##########    ########    ####    ###########       ############################
// ####################    #######            ######    ##########    ########    ####    ###########    ###############################
// ####################    #######            ######    ##########    ########    ####    ###########    ###############################
// ####################    #######    ####    ######    ##########    ########    ####    #####    ##    ###############################
// ####################         ##    ####    ######    ##########    ######        ###           ###         ##########################
// ####################         ##    ####    ######    ##########    ######        #####       #####         ##########################
// #####################################################################################################################################
// #####################################################################################################################################









/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/



