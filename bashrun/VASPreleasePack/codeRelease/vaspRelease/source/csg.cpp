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
 * File: csg.cpp
 *       Implementation for CSG operations
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

#include "csg.h"

// #####################################################################################################################################
// #####################################################################################################################################
// ######################       ########      ########       ###########################################################################
// ####################           ####          ####           #########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ###################    ###########    ##########    #################################################################################
// ###################    ############     ########    #################################################################################
// ###################    ##############      #####    ###      ########################################################################
// ###################    ##################    ###    ###      ########################################################################
// ###################    ###################    ##    #####    ########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ####################           ####          ####           #########################################################################
// ######################       ########      ########       ###########################################################################
// #####################################################################################################################################
// #####################################################################################################################################

/////////////////////////////////////////////////////////////////////////////////////
// Nothing is computed in the constructor.  Just set up the data.
/////////////////////////////////////////////////////////////////////////////////////
CSG::CSG( SurfaceObject * surf1, SurfaceObject * surf2 )
{
	
	// This is provided in runCSG()
	//op = CSG_op;

	runnable = true;
	lat1 = NULL;
	lat2 = NULL;

	if( surf1->numPoints == 0 || surf1->numTriangles == 0 ){
		printf("WARN! [CSG::CSG()           ]: Surface 1 has no geometry. \n");
		runnable = false;
	}
	if( surf2->numPoints == 0 || surf2->numTriangles == 0 ){
		printf("WARN! [CSG::CSG()           ]: Surface 2 has no geometry. \n");
		runnable = false;
	}

	if( runnable ){
		lat1 = new Lattice( surf1 );
		lat2 = new Lattice( surf2 );
	}

	
	edges = alloc_set(SP_MAP);
	edgePoints = alloc_set(SP_MAP);
}



/////////////////////////////////////////////////////////////////////////////////////
// destructor
/////////////////////////////////////////////////////////////////////////////////////
CSG::~CSG()
{
	int i = 0;
	int j = 0;

	// free the input lattices
	if( lat1 != NULL ){ delete(lat1); }
	if( lat2 != NULL ){ delete(lat2); }

	//This is a set of sets indicating the connection between cornerIndexes and other corner indexes
	for(i = 0; i<size_set(edges); i++){
		set_t thisSet = (set_t) mapsto_set(edges, edges[i]);
			for(j = 0; j<size_set(thisSet); j++){
				int * pt = (int *) mapsto_set(thisSet, thisSet[j]);
				delete[](pt);
			}
		free_set(thisSet);
	}
	free_set(edges);

	//This is a set of double*'s that are indexed by edges
	for(i = 0; i<size_set(edgePoints); i++){
		double * pt = (double *) mapsto_set(edgePoints, edgePoints[i]);
		delete[](pt);
	}
	free_set(edgePoints);
}



/////////////////////////////////////////////////////////////////////////////////////
// executive control function
/////////////////////////////////////////////////////////////////////////////////////
SurfaceObject * CSG::runCSG( int CSG_op, double resolution )
{
	SurfaceObject * result;
	res = resolution;
	op = CSG_op;

	// special case: one or both of the surfaces is empty.
	// we will try to return the right thing here.
	if (!runnable){
		result = degenerateProcessing();
		return result;
	}
	
	//get bounds
	lat1->setBounds(res);
	lat2->setBounds(res);

	//syncrhonize
	synchronize( lat1, lat2 );

	//now compute sign generation
	printf("########################## Generating Signed Field for Lattice 1 ##########################\n");
	lat1->insertTriangles();
	lat1->classifyCubes();
	lat1->fillNeighboringSets();
	lat1->populateFlags();
	
	printf("########################## Generating Signed Field for Lattice 2 ##########################\n");
	lat2->insertTriangles();
	lat2->classifyCubes();
	lat2->fillNeighboringSets();
	lat2->populateFlags();
	
	printf("############################### Now Computing CSG Operation ###############################\n");
	result = finalizeCSGsurface( );

	return result;
}











// ##########################################################################################################################
// ###     ###    ###     ##########     ###     ##  ###  #      #     ###     ##  ##########################################
// ##  ###  #  ##  #  ###  ########  ###  #  ###  #   ##  ###  ###  ##  #  ###  #  ##########################################
// ##  ######  #####  ###  ########  ###  #  ###  #    #  ###  ###  ##  #  ###  #  ##########################################
// ##  #######  ####  #############  ######  ###  #  #    ###  ###  ##  #  ###  #  ##########################################
// ##  #########  ##  #############  ######  ###  #  ##   ###  ###    ###  ###  #  ##########################################
// ##  ##########  #  ##   ########  ######  ###  #  ###  ###  ###  #  ##  ###  #  ##########################################
// ##  ##########  #  ###  ########  ###  #  ###  #  ###  ###  ###  ##  #  ###  #  ##########################################
// ##  ###  #  ##  #  ###  ########  ###  #  ###  #  ###  ###  ###  ##  #  ###  #  ##########################################
// ###     ###    ###     ##########     ###     ##  ###  ###  ###  ##  ##     ##     #######################################
// ##########################################################################################################################

// Functions for CSG capabilities.

/////////////////////////////////////////////////////////////////////////////////////
// Synchronize is a function that ensures that the lattices in both this lattice and
// "lat", supplied as input, have the same dimensions.  This ensures that point and
// cube indices one lattice map to the identical thing in the other lattice, so that
// we can make the interior/exterior assessments we need to make in both without 
// losing track of identity.  Call this function after setBounds() is called in both
// lattices, NEVER before.  op is the CSG operation to be used.
// This function also sets local variables (dimensions, resolution) for later use.
/////////////////////////////////////////////////////////////////////////////////////
void CSG::synchronize( Lattice * lat1, Lattice * lat2 )
{
	////////////////////////////////////////
	printf("STATE [synchronize()        ]: Synchronizing Lattice Bounds ... \n");
	////////////////////////////////////////

	if( lat1->res != lat2->res ){
		if( lat1->res < lat2->res ){ lat2->res = lat1->res; }
		else{ lat1->res = lat2->res; }
		printf("WARN  [synchronize()        ]: Resolutions don't match. Using the finer resolution (%f)\n", lat1->res );
	}

	///if this one is more extreme than the other one, set it's bounds to this.  otherwise do the reverse.
	if(lat1->xneg < lat2->xneg){ lat2->xneg = lat1->xneg; } else { lat1->xneg = lat2->xneg; }
	if(lat1->xpos > lat2->xpos){ lat2->xpos = lat1->xpos; } else { lat1->xpos = lat2->xpos; }
	if(lat1->yneg < lat2->yneg){ lat2->yneg = lat1->yneg; } else { lat1->yneg = lat2->yneg; }
	if(lat1->ypos > lat2->ypos){ lat2->ypos = lat1->ypos; } else { lat1->ypos = lat2->ypos; }
	if(lat1->zneg < lat2->zneg){ lat2->zneg = lat1->zneg; } else { lat1->zneg = lat2->zneg; }
	if(lat1->zpos > lat2->zpos){ lat2->zpos = lat1->zpos; } else { lat1->zpos = lat2->zpos; }

	if( lat1->xneg!=lat2->xneg || lat1->xpos!=lat2->xpos || lat1->yneg!=lat2->yneg ||
		lat1->ypos!=lat2->ypos || lat1->zneg!=lat2->zneg || lat1->zpos!=lat2->zpos)
	{
		printf("ERROR [synchronize()        ]: dimensions did not synchronize!\n"); exit(1);
	}

	//both these are already nudged in the negative direction, so we do not need a further nudge.
	///+1 for fractional cube round up (int conversion truncates)
	///= +1 total additions.
	int myxdim = ((int) ((lat1->xpos-lat1->xneg)/lat1->res)) + 1;       ///this is the total number of CUBES.  Not lines.
	int myydim = ((int) ((lat1->ypos-lat1->yneg)/lat1->res)) + 1;
	int myzdim = ((int) ((lat1->zpos-lat1->zneg)/lat1->res)) + 1;
	
	lat1->xdim = myxdim;	lat2->xdim = myxdim;
	lat1->ydim = myydim;	lat2->ydim = myydim;
	lat1->zdim = myzdim;	lat2->zdim = myzdim;

	///Reset the high dimensions based on the sizes of the cubes, since there is fractional roundup.
	lat1->xpos = lat1->xneg + (lat1->xdim*lat1->res);		lat2->xpos = lat2->xneg + (lat2->xdim*lat2->res);
	lat1->ypos = lat1->yneg + (lat1->ydim*lat1->res);		lat2->ypos = lat2->yneg + (lat2->ydim*lat2->res);
	lat1->zpos = lat1->zneg + (lat1->zdim*lat1->res);		lat2->zpos = lat2->zneg + (lat2->zdim*lat2->res);
        
	////////////////////////////////////////
	printf("INFO  [Lattice 1 Dimensions ]: xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", lat1->xneg, lat1->xpos, lat1->yneg, lat1->ypos, lat1->zneg, lat1->zpos);
	printf("INFO  [Lattice 1 Cube Dimens]: xdim: %i, ydim: %i, zdim: %i\n", lat1->xdim, lat1->ydim, lat1->zdim);
	printf("INFO  [Lattice 1 Resolution ]: %f\n", lat1->res);
	printf("INFO  [Lattice 1 Tot # Cubes]: %i\n", (lat1->xdim)*(lat1->ydim)*(lat1->zdim));

	printf("INFO  [Lattice 2 Dimensions ]: xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", lat2->xneg, lat2->xpos, lat2->yneg, lat2->ypos, lat2->zneg, lat2->zpos);
	printf("INFO  [Lattice 2 Cube Dimens]: xdim: %i, ydim: %i, zdim: %i\n", lat2->xdim, lat2->ydim, lat2->zdim);
	printf("INFO  [Lattice 2 Resolution ]: %f\n", lat2->res);
	printf("INFO  [Lattice 2 Tot # Cubes]: %i\n", (lat2->xdim)*(lat2->ydim)*(lat2->zdim));
	////////////////////////////////////////
	
	// set the local variables for later use.
	// resolution of the lattice
	res = lat1->res;

	// bounds of the lattice.
	xneg = lat1->xneg;	xpos = lat1->xpos;
	yneg = lat1->yneg;	ypos = lat1->ypos;
	zneg = lat1->zneg;	zpos = lat1->zpos;

	// dimensions of the lattice
	xdim = lat1->xdim;
	ydim = lat1->ydim;
	zdim = lat1->zdim;
}






/////////////////////////////////////////////////////////////////////////////////////
// This function generates the actual CSG surface.  It can only be called once sign
// generation on two synchronized lattices has been computed.  After that is done
// finalizeCSGsurface() finds all the nonempty cubes in the output surface (using
// helper function csgNonemptyCubeTest()) and then generates surface on them in the 
// same style as Lattice::finalizeSurface().  This data is then fed into the 
// SurfaceObject constructor, and the SurfaceObject is generated and returned.
/////////////////////////////////////////////////////////////////////////////////////
SurfaceObject * CSG::finalizeCSGsurface( )
{
	////////////////////////////////////////
	printf("STATE [finalizeCSGsurface() ]: Computing CSG output Surface ...");
	////////////////////////////////////////
	int i = 0;
	int j = 0;
	CubeTable * table = new CubeTable();
	
	// first, find the set of all nonempty cubes in both lattices.
	set_t testCubes = alloc_set(0);
	for(i = 0; i<size_set(lat1->nonempty); i++){
		testCubes = put_set(testCubes, lat1->nonempty[i]);
	}
	for(i = 0; i<size_set(lat2->nonempty); i++){
		testCubes = put_set(testCubes, lat2->nonempty[i]);
	}

	// Here we narrow it down to form the set of cubes that 
	// are nonempty in the output lattice
	set_t nonEmptyCubes = alloc_set(0);
	for(i = 0; i<size_set(testCubes); i++){
		int cubeIndex = testCubes[i];
		if( csgNonemptyCubeTest( cubeIndex ) ){
			nonEmptyCubes = put_set(nonEmptyCubes, cubeIndex );
		}
	}

	// the points are stored in edgePoints, so store the triangles here.
	set_t triList = alloc_set(SP_MAP);		//arrays of 3 ints, for triangle topology
	int * edgeHolder = new int[12];		//this array is used to map pointers to edge points.

	// Now, for each cube that will be nonempty in the final CSG object
	// we generate the surface.
	for(i = 0; i<size_set(nonEmptyCubes); i++){
		stdoutProgressBar(i, size_set(nonEmptyCubes));
		// clean up the edgeHolder
		for(j = 0; j<12; j++){ edgeHolder[j]=-1; }
		
		// first, get the cube's indices.	
		int cubeIndex = nonEmptyCubes[i];
		int cjkRemainder = cubeIndex % ((lat1->ydim)*(lat1->zdim));
		int ci = (int) ( ((double) (cubeIndex)) / ((double) ((lat1->ydim)*(lat1->zdim))  ) );
		int cj = (int) ( ((double) (cjkRemainder)) / ((double) (lat1->zdim)) );
		int ck = cjkRemainder % (lat1->zdim);
		
		// Name the corner indices						///get the corner status for this cube. 
		int c0 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+0) ;      int p0 = isInsideCSG( c0 );
		int c1 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+1) ;      int p1 = isInsideCSG( c1 );
		int c2 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+0) ;      int p2 = isInsideCSG( c2 );
		int c3 = (ci+0)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+1) ;      int p3 = isInsideCSG( c3 );
		int c4 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+0) ;      int p4 = isInsideCSG( c4 );
		int c5 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+0)*(zdim+1) + (ck+1) ;      int p5 = isInsideCSG( c5 );
		int c6 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+0) ;      int p6 = isInsideCSG( c6 );
		int c7 = (ci+1)*((ydim+1)*(zdim+1)) + (cj+1)*(zdim+1) + (ck+1) ;      int p7 = isInsideCSG( c7 );

		//If the corners have different interior/exterior settings, and an edge point has 
		//not been generated in the past, then store an edge point.
		if(p0!=p1){ edgeHolder[ 0] = getEdgePointCSG(c0,c1); if( edgeHolder[ 0] == -1){edgeHolder[ 0] = storeEdgePointCSG(c0, c1, computeEdgeCrossingPointCSG( c0, c1 ) ); } }
		if(p0!=p2){ edgeHolder[ 1] = getEdgePointCSG(c0,c2); if( edgeHolder[ 1] == -1){edgeHolder[ 1] = storeEdgePointCSG(c0, c2, computeEdgeCrossingPointCSG( c0, c2 ) ); } }
		if(p0!=p4){ edgeHolder[ 4] = getEdgePointCSG(c0,c4); if( edgeHolder[ 4] == -1){edgeHolder[ 4] = storeEdgePointCSG(c0, c4, computeEdgeCrossingPointCSG( c0, c4 ) ); } }
		if(p4!=p5){ edgeHolder[ 8] = getEdgePointCSG(c4,c5); if( edgeHolder[ 8] == -1){edgeHolder[ 8] = storeEdgePointCSG(c4, c5, computeEdgeCrossingPointCSG( c4, c5 ) ); } }
		if(p4!=p6){ edgeHolder[ 9] = getEdgePointCSG(c4,c6); if( edgeHolder[ 9] == -1){edgeHolder[ 9] = storeEdgePointCSG(c4, c6, computeEdgeCrossingPointCSG( c4, c6 ) ); } }
		if(p2!=p3){ edgeHolder[ 3] = getEdgePointCSG(c2,c3); if( edgeHolder[ 3] == -1){edgeHolder[ 3] = storeEdgePointCSG(c2, c3, computeEdgeCrossingPointCSG( c2, c3 ) ); } }
		if(p2!=p6){ edgeHolder[ 6] = getEdgePointCSG(c2,c6); if( edgeHolder[ 6] == -1){edgeHolder[ 6] = storeEdgePointCSG(c2, c6, computeEdgeCrossingPointCSG( c2, c6 ) ); } }
		if(p1!=p3){ edgeHolder[ 2] = getEdgePointCSG(c1,c3); if( edgeHolder[ 2] == -1){edgeHolder[ 2] = storeEdgePointCSG(c1, c3, computeEdgeCrossingPointCSG( c1, c3 ) ); } }
		if(p1!=p5){ edgeHolder[ 5] = getEdgePointCSG(c1,c5); if( edgeHolder[ 5] == -1){edgeHolder[ 5] = storeEdgePointCSG(c1, c5, computeEdgeCrossingPointCSG( c1, c5 ) ); } }
		if(p6!=p7){ edgeHolder[11] = getEdgePointCSG(c6,c7); if( edgeHolder[11] == -1){edgeHolder[11] = storeEdgePointCSG(c6, c7, computeEdgeCrossingPointCSG( c6, c7 ) ); } }
		if(p3!=p7){ edgeHolder[ 7] = getEdgePointCSG(c3,c7); if( edgeHolder[ 7] == -1){edgeHolder[ 7] = storeEdgePointCSG(c3, c7, computeEdgeCrossingPointCSG( c3, c7 ) ); } }
		if(p5!=p7){ edgeHolder[10] = getEdgePointCSG(c5,c7); if( edgeHolder[10] == -1){edgeHolder[10] = storeEdgePointCSG(c5, c7, computeEdgeCrossingPointCSG( c5, c7 ) ); } }

		///figure out what case we in the cubetable we will use.
		int cubeTableIndex = p0 + 2*p1 + 4*p2 + 8*p3 + 16*p4 + 32*p5 + 64*p6 + 128*p7;
		int numPolysToAdd = table->sizes[cubeTableIndex];
		int * polyList = table->polys[cubeTableIndex];

		for(j = 0; j<numPolysToAdd; j++){
			int * newTriangle = new int[3];
			///reversed order here to reverse norms - original polyList is reversed.
			newTriangle[0] = edgeHolder[ polyList[3*j+0] ];
			newTriangle[1] = edgeHolder[ polyList[3*j+2] ];
			newTriangle[2] = edgeHolder[ polyList[3*j+1] ];

			///store the triangle topology.
			triList = associate_set(triList, size_set(triList), newTriangle);
		}
	}

//	////////////////////////////////////////
	printf("INFO  [# Geo. Tested Points ]: %i \n", lat1->numGeometricChecks + lat1->numGeometricChecks);
	printf("INFO  [Total Pts processed  ]: %i \n", size_set(lat1->points) + size_set(lat1->points));
	printf("INFO  [finalizeSurface()    ]: Generating Surface... \n");
//	////////////////////////////////////////

	SurfaceObject * result = new SurfaceObject(edgePoints, triList);

	////////////////////////////////////////
	printf("INFO  [finalizeSurface()    ]: Done generating surface. \n");
	////////////////////////////////////////

	free_set(testCubes);
	free_set(nonEmptyCubes);
	delete(table);
	delete[](edgeHolder);
	for(i = 0; i<size_set(triList); i++){
		int * tri = (int *) mapsto_set(triList, i);
		delete[](tri);
	}
	free_set(triList);
	
	return result;
}





/////////////////////////////////////////////////////////////////////////////////////
// This function processes the case where either surface actually contains no geometry.
// In these cases, we may return empty, or we may return one of the surfaces, based on
// the operation and on the inputs.  
/////////////////////////////////////////////////////////////////////////////////////
SurfaceObject * CSG::degenerateProcessing( )
{
	SurfaceObject * result;
	SurfaceObject * s1 = NULL;
	SurfaceObject * s2 = NULL;

	if( lat1 != NULL ){ s1 = lat1->surf; }
	if( lat2 != NULL ){ s2 = lat2->surf; }
	
	if( op == CSG_INTERSECTION ){
		result = new SurfaceObject();
	}	
	if( op == CSG_DIFFERENCE ){
		if( s1 == NULL && s2 == NULL ){
			result = new SurfaceObject();
		}
		if( s1 == NULL && s2 != NULL ){
			result = new SurfaceObject();
		}
		if( s1 != NULL && s2 == NULL ){
			result = s1->copy();
		}
	}
	if( op == CSG_UNION ){
		if( s1 == NULL && s2 == NULL ){
			result = new SurfaceObject();
		}
		if( s1 == NULL && s2 != NULL ){
			result = s2->copy();
		}
		if( s1 != NULL && s2 == NULL ){
			result = s1->copy();
		}
	}
	
	return result;
}









// ####################################################################################################################################################
// ###     ###    ###     #####  ###  #      #  #####     ##      #     #####      #  ##  #  ###  ##    ##      #    ##    ##  ###  ##    #############
// ##  ###  #  ##  #  ###  ####  ###  #  #####  #####  ##  #  #####  ##  ####  #####  ##  #   ##  #  ##  ###  ####  ##  ##  #   ##  #  ##  ############
// ##  ######  #####  ###  ####  ###  #  #####  #####  ##  #  #####  ##  ####  #####  ##  #    #  #  #######  ####  ##  ##  #    #  #  ################
// ##  #######  ####  #########       #    ###  #####  ##  #    ###  ##  ####     ##  ##  #  #    #  #######  ####  ##  ##  #  #    ##  ###############
// ##  #########  ##  #########  ###  #  #####  #####     ##  #####    ######  #####  ##  #  ##   #  #######  ####  ##  ##  #  ##   ####  #############
// ##  ##########  #  ##   ####  ###  #  #####  #####  #####  #####  #  #####  #####  ##  #  ###  #  #######  ####  ##  ##  #  ###  #####  ############
// ##  ##########  #  ###  ####  ###  #  #####  #####  #####  #####  ##  ####  #####  ##  #  ###  #  #######  ####  ##  ##  #  ###  #####  ############
// ##  ###  #  ##  #  ###  ####  ###  #  #####  #####  #####  #####  ##  ####  #####  ##  #  ###  #  ##  ###  ####  ##  ##  #  ###  #  ##  ############
// ###     ###    ###     #####  ###  #      #      #  #####      #  ##  ####  ######    ##  ###  ##    ####  ###    ##    ##  ###  ##    #############
// ####################################################################################################################################################

// Functions that support CSG capabilities

/////////////////////////////////////////////////////////////////////////////////////
// This function determins if a point is inside or outside based on the
// desired CSG operation, and the interior/exterior state of the cubes
// in each Lattice object.
// This returns only 1 or 0, inside or outside.
/////////////////////////////////////////////////////////////////////////////////////
int CSG::isInsideCSG( int pointIndex )
{
	int result = -1;

	if( op == CSG_INTERSECTION ){
		//Perform rapid caching tests first.
		int quickTest1 = lat1->getPointStatus(pointIndex);
		if ( quickTest1 == LATTICE_EXTERIOR ){ return 0; }
		int quickTest2 = lat2->getPointStatus(pointIndex);
		if ( quickTest2 == LATTICE_EXTERIOR ){ return 0; }
		if( quickTest1 == LATTICE_INTERIOR && quickTest2 == LATTICE_INTERIOR ){ return 1; }

		//Caching did not have enough to answer the question, so perform a full test.
		if ( quickTest1 != LATTICE_INTERIOR && quickTest1 != LATTICE_EXTERIOR ){
			int test1 = lat1->isInside( pointIndex );
			if( test1==0 ){ return 0; }
			else{ quickTest1 = LATTICE_INTERIOR; }
		}
		if ( quickTest2 != LATTICE_INTERIOR && quickTest2 != LATTICE_EXTERIOR ){
			int test2 = lat2->isInside( pointIndex );
			if( test2==0 ){ return 0; }
			else{ quickTest2 = LATTICE_INTERIOR; }
		}
		
		if( quickTest1 == LATTICE_INTERIOR && quickTest2 == LATTICE_INTERIOR ){ result = 1; }
		else{ result = 0; }
	}
	
	if( op == CSG_UNION ){
		//Perform rapid caching tests first.
		int quickTest1 = lat1->getPointStatus(pointIndex);
		if ( quickTest1 == LATTICE_INTERIOR ){ return 1; }
		int quickTest2 = lat2->getPointStatus(pointIndex);
		if ( quickTest2 == LATTICE_INTERIOR ){ return 1; }
		if( quickTest1 == LATTICE_EXTERIOR && quickTest2 == LATTICE_EXTERIOR ){ return 0; }

		//Caching did not have enough to answer the question, so perform a full test.
		if ( quickTest1 != LATTICE_INTERIOR && quickTest1 != LATTICE_EXTERIOR ){
			int test1 = lat1->isInside( pointIndex );
			if( test1==1 ){ return 1; }
			else{ quickTest1 = LATTICE_EXTERIOR; }
		}
		if ( quickTest2 != LATTICE_INTERIOR && quickTest2 != LATTICE_EXTERIOR ){
			int test2 = lat2->isInside( pointIndex );
			if( test2==1 ){ return 1; }
			else{ quickTest2 = LATTICE_EXTERIOR; }
		}

		if( quickTest1 == LATTICE_EXTERIOR && quickTest2 == LATTICE_EXTERIOR ){ result = 0; }
		else{ result = 1; }
	}
	
	if( op == CSG_DIFFERENCE ){
		//Perform rapid caching tests first.
		int quickTest1 = lat1->getPointStatus(pointIndex);
		if ( quickTest1 == LATTICE_EXTERIOR ){ return 0; }
		int quickTest2 = lat2->getPointStatus(pointIndex);
		if ( quickTest2 == LATTICE_INTERIOR ){ return 0; }
		if( quickTest1 == LATTICE_INTERIOR && quickTest2 == LATTICE_EXTERIOR ){ return 1; }

		//Caching did not have enough to answer the question, so perform a full test.
		if ( quickTest1 != LATTICE_INTERIOR && quickTest1 != LATTICE_EXTERIOR ){
			int test1 = lat1->isInside( pointIndex );
			if( test1==0 ){ return 0; }
			else{ quickTest1 = LATTICE_INTERIOR; }
		}
		if ( quickTest2 != LATTICE_INTERIOR && quickTest2 != LATTICE_EXTERIOR ){
			int test2 = lat2->isInside( pointIndex );
			if( test2==1 ){ return 0; }
			else{ quickTest2 = LATTICE_EXTERIOR; }
		}

		if( quickTest1 == LATTICE_INTERIOR && quickTest2 == LATTICE_EXTERIOR ){ result = 1; }
		else{ result = 0; }
	}

	if( result == -1 ){
		printf("ERROR [isInsideCSG()  ]: Function did not decide CSG status!\n");
		exit(1);
	}
	return result;
}





/////////////////////////////////////////////////////////////////////////////////////
// This function determines if a given cube is nonempty in the output
// surface, given the two input surfaces and the desired CSG operation.
// Returns true if this should be nonempty in the output,
// returns false if this should be empty in the output.
/////////////////////////////////////////////////////////////////////////////////////
bool CSG::csgNonemptyCubeTest( int cubeIndex )
{
	bool inLat1 = contains_set( lat1->nonempty, cubeIndex );
	bool inLat2 = contains_set( lat2->nonempty, cubeIndex );

	//if it is empty in both, then it's empty in the output, for all operations. (should never happen)
	if( !inLat1 && !inLat2 ){ 
		printf("ERROR [csgNonemptyCubeTest()]: This should never happen!\n");
		return false; 
	}

	//if it is nonempty in both, then it must be analyzed for geometry.  so it is "nonempty"
	if(  inLat1 &&  inLat2 ){ return true;  }

	///from here on, one side is nonempty, one side is empty.
	bool result;

	///process the case of intersection
	if( op == CSG_INTERSECTION ){
		// find the state of the empty cube.
		int flagVal;
		if( !inLat1 ){ flagVal = lat1->getEmptyCubeStatus(cubeIndex); }
		if( !inLat2 ){ flagVal = lat2->getEmptyCubeStatus(cubeIndex); }
		if( flagVal != LATTICE_INTERIOR && flagVal != LATTICE_EXTERIOR ){
			printf("ERROR [csgNonemptyCubeTest()]: found an unflagged empty cube! flags must be found first!\n");
			exit(1);
		}

		//for INTERSECTION, if this geometry is in the interior of the other, then its relevant
		if( flagVal == LATTICE_INTERIOR ){ return true; }
		//otherwise it is in the exterior, in which case it is not in the intersection, so eliminate
		else{ result = false; }
	}

	///process the case of union
	if( op == CSG_UNION ){
		// find the state of the empty cube.
		int flagVal;
		if( !inLat1 ){ flagVal = lat1->getEmptyCubeStatus(cubeIndex); }
		if( !inLat2 ){ flagVal = lat2->getEmptyCubeStatus(cubeIndex); }
		if( flagVal != LATTICE_INTERIOR && flagVal != LATTICE_EXTERIOR ){
			printf("ERROR [csgNonemptyCubeTest()]: found an unflagged empty cube! flags must be found first!\n");
			exit(1);
		}

		//for UNION, if this geometry is in the exterior of the other, then its relevant
		if( flagVal == LATTICE_EXTERIOR ){ return true; }
		//otherwise it is in the interior, in which case it is does not define the UNION border
		else{ result = false; }
	}

	///process the case of difference
	//note the second lattice is always subtracting away from the first.
	if( op == CSG_DIFFERENCE ){
		// find the state of the empty cube.
		int flagVal;
		if( !inLat1 ){ 
			flagVal = lat1->getEmptyCubeStatus(cubeIndex); 
			if( flagVal != LATTICE_INTERIOR && flagVal != LATTICE_EXTERIOR ){
				printf("ERROR [csgNonemptyCubeTest()]: found an unflagged empty cube! flags must be found first!\n");
				exit(1);
			}
			if( flagVal == LATTICE_INTERIOR ){ return true; }
			else{ result = false; }
		}
		if( !inLat2 ){ 
			flagVal = lat2->getEmptyCubeStatus(cubeIndex);
			if( flagVal != LATTICE_INTERIOR && flagVal != LATTICE_EXTERIOR ){
				printf("ERROR [csgNonemptyCubeTest()]: found an unflagged empty cube! flags must be found first!\n");
				exit(1);
			}
			if( flagVal == LATTICE_INTERIOR ){ return false; }
			else{ result = true; }
		}
	}
	
	return result;

}




/////////////////////////////////////////////////////////////////////////////////////
// This function searches the nested set structure "edges" to determine if it contains
// an index corresponding to c0 and c1 for a geometric point in edgePoints.  
// If if a point exists, the function returns the index, otherwise it returns -1.
// Note this is basically identical to the fucntion in Lattice
/////////////////////////////////////////////////////////////////////////////////////
int CSG::getEdgePointCSG(int p0, int p1)
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
// This function inserts an edge into "edges", for indices c0 and c1, that corresponds
// to point "pt".  "pt" itself is stored in edgePoints, using the index stored in edges.
// Returns the index into edgePoints for this edge POint.
// Note this is basically identical to the fucntion in Lattice
/////////////////////////////////////////////////////////////////////////////////////
int CSG::storeEdgePointCSG(int c0, int c1, double * pt )
{
	///a must always be lower, so we can save space.
	int a = c0;
	int b = c1;
	if(a > b){
		a = c1;
		b = c0;
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
// This function computes the point at which the CSG boundary intersects the edge of
// the cube.  This is achieved by first getting the interior/exterior status of the
// points in each Lattice, then making a decision based on the appropriate intersection
// point.
/////////////////////////////////////////////////////////////////////////////////////
double * CSG::computeEdgeCrossingPointCSG( int p0, int p1 )
{
	// first get the inside/outside state.  These calls were all made before,
	// so they will not be slow or an excessive use of CPU.
	int p0inLat1 = lat1->isInside( p0 );
	int p0inLat2 = lat2->isInside( p0 );
	int p1inLat1 = lat1->isInside( p1 );
	int p1inLat2 = lat2->isInside( p1 );

	// first, if there is only one intersection, return the intersection.  Trivial.
	if( !(p0inLat1==p1inLat1) && (p0inLat2==p1inLat2) ){ return lat1->getAverageIntersectionPt(p0, p1); }
	if( (p0inLat1==p1inLat1) && !(p0inLat2==p1inLat2) ){ return lat2->getAverageIntersectionPt(p0, p1); }

	// If we got to this point, that means that both intersect.
	// One of the points is inside the CSG surface, one is outside.
	// we can find out, because it is cached now.
	int p0inCSG = isInsideCSG( p0 );
	int insidePt;
	if( p0inCSG ){ insidePt = p0; }
	else{ insidePt = p1; }

	// since both segments intersect the CSG surface find out 
	// the points at which they intersect the surface.
	double * pointa = lat1->getAverageIntersectionPt(p0, p1);	//the intersect on Lat1
	double * pointb = lat2->getAverageIntersectionPt(p0, p1);	//the intersect on Lat2

	// get the coordinates of the inside point.
	double * inPt = lat1->getPoint(insidePt);
	double dista = vectorSize(inPt, pointa);
	double distb = vectorSize(inPt, pointb);

	double * output;	//this just points to the answer, doesnt store it.

	// handle each case:
	// In the union, the point further from the interior point is 
	// the border of the union
	if( op == CSG_UNION ){
		if( dista > distb ){ output = pointa; }
		else{ output = pointb; }
	}
	
	// In the Intersection, the point closer to the interior point
	// is the border of the intersection
	if( op == CSG_INTERSECTION ){
		if( dista < distb ){ output = pointa; }
		else{ output = pointb; }
	}

	// In the Difference, the point closer to the interior point
	// is the border of the difference
	if( op == CSG_DIFFERENCE ){
		if( dista < distb ){ output = pointa; }
		else{ output = pointb; }
	}
	
	//new the result 
	double * result = new double[3];
	result[0] = output[0];
	result[1] = output[1];
	result[2] = output[2];
	
	delete[](pointa);
	delete[](pointb);
	delete[](inPt);
	
	return result;

}









// #####################################################################################################################################
// #####################################################################################################################################
// ######################       ########      ########       ###########################################################################
// ####################           ####          ####           #########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ###################    ###########    ##########    #################################################################################
// ###################    ############     ########    #################################################################################
// ###################    ##############      #####    ###      ########################################################################
// ###################    ##################    ###    ###      ########################################################################
// ###################    ###################    ##    #####    ########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ####################           ####          ####           #########################################################################
// ######################       ########      ########       ###########################################################################
// #####################################################################################################################################
// #####################################################################################################################################




/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/



