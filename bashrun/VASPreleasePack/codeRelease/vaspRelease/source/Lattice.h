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
 * File: Lattice.h
 *       interface for the Lattice structure, for managing marching cubes
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VOLUMETRIC ANALYSIS OF SURFACE PROPERTIES (VASP)                                                           //
// BRIAN Y CHEN, HONIG LAB, 2010                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// UnifiedMarchingCubes.h: interface for simplification instructions
//
//////////////////////////////////////////////////////////////////////
#ifndef UNIFIED_MARCHING_CUBES
#define UNIFIED_MARCHING_CUBES

#define LATTICE_INTERIOR 100
#define LATTICE_EXTERIOR 200
#define LATTICE_NONEMPTY 300
#define LATTICE_NOT_SET  400
#define LATTICE_RADIUS_TO_TEST 2

//make this odd, so that voting does not end in tie.
#define NUMBER_OF_CONTINGENCY_TESTS 51
#define NUMBER_OF_CROSS_SET_TESTS 11

#include "StdAfx.h"

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
// For sets specifying cube neighbors, they are specified [high-x][low-x][high-y][low-y][high-z][low-z]
// In the case where the cube has no neighbor in a certain position, NULL is provided in that space.
// Thus, sets specifying cube neighbors always have exactly six elements, but may have as few as three
// non-NULL members.

// Design Note:
// THis class is designed with the idea that there exist a set of points, a set of edges, and a set
// of cubes, all independently represented and correctly marked as adjacent.  However, since this
// is a uniform grid, I can procedurally generate most of this information, and thus only keep a 
// tiny amount to actually remember.  Thus, there were originally "LatticePOints" "LatticeEdges" 
// and "LatticeCubes" that are now unnecessary.

class Lattice
{

	public:

	// constructor and destructor
	Lattice( SurfaceObject * s );
	virtual ~Lattice();

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

	// Primary Control methods.  To be run in this order.
       
	// setBounds() sets the boundaries of the Lattice, given the surface (surf)
	// stored in the Lattice ONLY.  For CSG operations with multiple lattices,
	// The lattices must be synchronized using an external function (FIXME)
	void setBounds(double r);
	
	// Scan conversion: 
	// This function inserts all triangles from the surfaceObject into cubes.
	// The cubes are virtual.  They are conceptually indexed in a grid, with indices
	//    i*(ydim*zdim) + j*(zdim) + k;
	void insertTriangles();
	
	// Cube sign generation step 1
	// This function finds all the connected components of empty and nonempty cubes.
	// Lists of indices of empty cubes go into emptyGroups
	// Lists of indices of nonempty cubes go into nonemptyGroups
	// Allocates the flags array, which could not be allocated without the number of
	// empty groups
	void classifyCubes();

	// Cube sign generation step 2
	// This function fills the neighborMatrix, with segment-based testing on
	// the nonempty cubesets.
	void fillNeighboringSets();

	// Cube sign generation step 3 
	// This function runs segment-based tests on adjacent empty groups, based on the information
	// in the neighborMatrix.  The neighbormatrix specifies specific cubes that are adjacent to
	// two different groups of empty connected cubes.  Segment tests are run in the vicinity of
	// those specific cubes, to ensure that the segments tested are as short as possible.
	// Furthermore, using breath first search, tests are only run as needed to propagate inside-
	// outside information from the exterior (which is, by definition, exterior) too all other
	// empty groups.  As a result except for some special cases, all empty groups should be
	// properly flagged in using this function.  flagProximalSets covers all remaining cases.
	void populateFlags();
	
	// Cube sign generation step 4 (contingency - not always called)
	// This function is called only when fillNeighboringSets() is unable to assign a flag to all
	// connected components of empty cubes based on neighborhood adjacency.  This function uses the
	// proximityTable (which is filled by fillNeighboringSets()) to then run cube-to-cube segment 
	// tests and assign interior/exterior values.
	void flagRemainingSets();

	// This function aggregates all the point values, either from caching, cube-flagging, or from
	// segment testing.  Then it calls getAverageIntersectionPt() and computes intersection
	// points for every edge along two different point values, on a nonempty cube, and stores it
	// so that it never needs to be computed again.  Then, on each nonempty cube, it computes
	// triangles.  Then it initializes the output surface, and returns.
	SurfaceObject * finalizeSurface();



	
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

	// Heavy Geometric Support Functions//////////////////////////////////////////////////////

	// Wrapper function for interior/exterior testing.
	// First calls getPointStatus, to check information based on empty cubes.
	// This should get the vast majority of points.
	// Then calls testShortRay, to do a segment test if the point is totally internal
	// Note that this MUST return an answer, and cannot return anything except 
	// 1, for inside OR 0 for outside
	int isInside( int pointIndex );

	// shortRay Test:
	// This function finds nearby empty cubes, generates a segment to those cubes, and 
	// then runs countIntersectionPts to count the number of intersections between the point
	// and the cube.  The result is a poll of the number of local cubes.
	// the result is returned as a value from 0 to 1, where 0.0 is external, and 1.0 is internal,
	// and any fraction between 0 and 1 is fraction of votes preferring internal relative to external
	double testShortRay(int pointIndex);

	// Segment testing:
	// This function takes a segment and tests how many times it nontrivially intersects
	// the input surface.  This eliminates double-intersections where it strikes exactly
	// at the edge between two triangles with an epsilon test - thus, the points must
	// be apart from each other.  Note that it does not catch the (very very unlikely) chance
	// that the segment intersects the segment at which two triangles intersect, where the
	// triangles do not pass through the segment - i.e. hitting the top of a mountain.
	int countIntersectionpts(double * a, double * b);
	
	// Segment testing:
	// This function counts the sections along a segment between the centers of two
	// cubes (e.g. cubetoCubeSegmentTest( testCube, extCube ) ).  A wrapper for the
	// function offers the functionality where members of two groups are specified
	// as neighbors to a specific cube, and all ordered pairs of cubes, one from one
	// group (e.g. currentGroup) and one from the other (e.g. adjacentGroup) are tested
	// using cubetoCubeSegmentTest( testCube, extCube ).
	int cubetoCubeSegmentTest( int testCube, int extCube );
	int cubetoCubeSegmentTest( int currentGroup, int adjacentGroup, int connectingCube );

	// Segment testing:
	// This function gets the intersections along a single edge in the lattice.
	// the input are the indices of the two points along the cube to be tested
	// Note that even thouse the points are adjacent to up to four cubes, anything
	// that intersects the edge is in all four cubes.  So we get all four, and
	// eliminate anything that is not in all four.  Then we test segment-tri intersections
	// and get the position of the average intersection.
	double * getAverageIntersectionPt(int pt1, int pt2);


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

	// Helper Functions /////////////////////////////////////////////////////////////////

	// Helper Function - Given the index of a cube, returns the indices of
	// cubes adjacent to this one, or returns NULL if it is out of range.
	// The cubes are virtual.  They are conceptually indexed in a grid, with indices
	//    i*(ydim*zdim) + j*(zdim) + k;
	// will not include an adjacent cube in set_t "excluded", if this is non-null;
	// output is added to outputSet, which MUST NOT BE NULL.
	// if outputSet is NULL, allocate a new one.
	// Here adjacency applies ONLY TO CUBES THAT SHARE A FACE.
	set_t getAdjacentCubes( int index, set_t excluded, set_t outputSet );

	// Helper function - this gives the set of cubes that shares a face or a
	// point with the cube provided (by it's index), except for any excluded
	// cube.  Results are put into outputSet, unless outputSet is NULL, in
	// which case a new set is allocated and returned.
	// Here NEIGHBOR applies CUBES THAT SHARE A FACE OR CORNER.
	set_t getNeighborCubes( int index, set_t excluded, set_t outputSet );

	// Helper Function - Union find on cube adjacency.
	// This performs a rapid union find on the cubes that are indexed in set_t cubes.
	// The assumed ajacency rules are applied, except in the case of being adjacent
	// to something in excluded - those cannot be part of any union.
	// Returns a set of sets, which define the equivalence classes.
	set_t unionFind( set_t excluded );

	// Helper Function - gets the cubes that this segment intersects.  exact.
	// note that the segment will *almost* always start with point a on a cube corner.
	// THis is called recursively until it gets to the end.
	// IMPORTANT: this is only really meant to be called where the first point is on the edge of a cube.
	// not in the middle - there are nudging values that could nudge it out of the cube otherwise.
	// IMPORTANT: if the first point v1 starts out of bounds, then this will return nothing, even if the segment
	// eventually comes in bounds.
	set_t getSegmentCubes(double * v1, double * v2);

	// Helper Function - Finds out if the point in question is inside or outside based
	// on assigned cubes.  Beginning with the index of the point, which are indexed as follows,
	//    i*((ydim+1)*(zdim+1)) + j*(zdim+1) + k;
	// we then find the index of the cubes it is adjacent to, and check their int/ext status.
	// returns LATTICE_INTERIOR if the point is inside, LATTICE_EXTERIOR of the point is outside
	// and LATTICE_NONEMPTY if the point is surrounded by nonempties.
	int getPointStatus(int pointIndex);

	// Helper function finds the cubes adjacent to a point, given the point index.
	// This basically handles the indexing.  Radius 0 is nothing, radius 1 is all 8 cubes
	// that have this point as a corner, radius 1 is all cubes around those cubes, etc.
	set_t getCubesAdjacentToPoint(int pointIndex, int radius);

	// Helper function gets the double vector position of a point using it's index.
	double * getPoint(int index);
	
	// Helper function gets the flag assigned to a specific cube
	int getEmptyCubeStatus(int cubeIndex);

	// Helper functions get and store edgePoints into edges
	// Returns the index into edgePoints for this edge POint.
	int storeEdgePoint(int p0, int p1, double * pt);
	// returns NULL if point pair is not stored.
	int getEdgePoint(int p0, int p1);

	// helper function returns a set containing the indices of the empty groups that
	// have cubes adjacent to the input cube
	set_t getNeighboringEmptyGroups( int cubeIndex);


	// ##########################################################################################################################
	// ##  ###  #### ####     ##    #### ####     ##  #####      ##    ##########################################################
	// ##  ###  ###   ###  ##  ##  ####   ###  ##  #  #####  #####  ##  #########################################################
	// ##  ###  ##     ##  ##  ##  ###     ##  ##  #  #####  #####  #############################################################
	// ##  ###  #   #   #  ##  ##  ##   #   #     ##  #####    ####  ############################################################
	// ##  ###  #  ###  #    ####  ##  ###  #  ##  #  #####  ########  ##########################################################
	// ##  ###  #       #  #  ###  ##       #  ##  #  #####  #########  #########################################################
	// ###  #  ##  ###  #  ##  ##  ##  ###  #  ##  #  #####  #########  #########################################################
	// ####   ###  ###  #  ##  ##  ##  ###  #  ##  #  #####  #####  ##  #########################################################
	// ##### ####  ###  #  ##  #    #  ###  #     ##      #      ##    ##########################################################
	// ##########################################################################################################################
	
	//Variables/////////////////////////////////////////////////////////////////

	double res;							//resolution of the lattice
	double xneg, xpos, yneg, ypos, zneg, zpos;	//bounds of the lattice.
	int xdim, ydim, zdim;					//number of cubes on each side.
	
	SurfaceObject * surf;					// Surface that this lattice describes.

	set_t nonempty;						// Set pointing to cubes that contain geometry.
										// This array is indexed as described in the	
										// insertTriangles() comment above.

	set_t nonemptyGroups;					// A set of sets listing indices of connected 
										// components of adjacent non-empty cubes.

	set_t emptyGroups;						// A set of sets listing indices of connected 
										// components of adjacent empty cubes.

	set_t neighborMatrix;					// A set of sets of ints.  Top level indices correspond 
										// to EMPTY connected components and link to empty sets, 
										// whose (second-level) indices correspond to other empty
										// sets, and link to an int[1], which stores the index
										// of the cube that is adjacent to cubes of both empty
										// sets.

	int * flags;							// int array stores flaggs (i.e. LATTICE_INTERIOR, etc)
										// for the empty cube groups.;

	set_t edges;							// This set contains sets that map to a point in space.
	set_t edgePoints;						// The top level set is indexed by points, and the
										// lower level set is also indexed by points, so you get
										// a way to index point-to-point connections - i.e. edges.
										// The lower level set contains int* values (int[1]) that index
										// into edgePOints, where you can get the actual double vector

	set_t points;							// This set contains pointers to point_inside or
	int * point_inside;						// point_outside, indexed on point indices.  this simply
	int * point_outside;					// states if a given point is inside or outside, and is only
										// applied to points that are not already set by cube
										// flagging.

	int numGeometricChecks;					// This diagnostic variable stores the number of times we had to
										// determine interior/exterior checking via segment testing, and
										// not based on cube flagging.
};




#endif


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







