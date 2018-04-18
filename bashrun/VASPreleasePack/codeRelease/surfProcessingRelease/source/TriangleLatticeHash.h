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
 * File: TriangleLatticeHash.h
 *       This class generates a rapid retrieval data structure for triangles in connolly surfaces,
 *       binding each triangle to the cubes in the lattice that it intersects with.  This permits points
 *       to determine rapidly if they are in the interior or exterior of the surface, and also to
 *       determine how far they are from the surface.
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

#ifndef TRIANGLE_LATTICE_HASH_H
#define TRIANGLE_LATTICE_HASH_H

#include "StdAfx.h"

///this  is how much extra space should be stored in the lattice outside of the max size of the surface.
#define LATTICE_PADDING 5

///we will check inside/outside status based on segment tests to nearby cubes.
///LOCAL_CUBE_CHECK_MULTIPLIER determines the radius of nearby cubes to check.
#define LOCAL_CUBE_CHECK_MULTIPLIER 1.5

///We will check at most this many cubes, for local cube checking.
#define LOCAL_CUBE_CHECK_MAX_CUBES 5



/*
#########################################################################################
#   #   ##   ##  ###  ## ###   # ####  ##   #   #   ##  ##   # ## ##  ###  ## ## ########
## ## ## ## ## ## # ## # ### ### ### ## ## ### ### ## ## # ### ## # ## # ## # ## ########
## ## ## ## ## ## # #### ###  ## ### ## ## ### ### ## ####  ##    # ## ## ###    ########
## ##   ### ##    # #  # ### ### ###    ## ### ### ## #### ### ## #    ### ## ## ########
## ## ## ## ## ## # ## # ### ### ### ## ## ### ### ## ## # ### ## # ## # ## # ## ########
## ## ## #   # ## ##  ##   #   #   # ## ## ### ##   ##  ##   # ## # ## ##  ## ## ########
#########################################################################################
this code uses the following cube layout, with 3D axes to the right.
      7---6            ^
      |\  |\       \   |
      | 4---5       +y |
      | | | |        \ z+
      3-|-2 |         \|  
       \|  \|          0--x+-->    
        0---1
*/

class LatticeObj;
class SurfaceObject;

///This lightweight class keeps a nonredundant list of 3d points.
///It tests for nonredundancy by linear search, so it cannot
///    handle very many points without becoming slow.
class SmallPointHash
{
	public:

	SmallPointHash();
	virtual ~SmallPointHash();

	set_t points;
	
	void addPoint(double * pt);
	void addPoint(double x, double y, double z);
	bool testPoint(double x, double y, double z);
	int size();
	double * getPoint(int index);
};




class TriangleLatticeHash  
{
	public:

	TriangleLatticeHash(SurfaceObject * sobj, double resolution, char * n);
	virtual ~TriangleLatticeHash();

	///returns the related cube for x cubes along the x dir, y cubes along the y dir, z cubes along the z dir
	/// if the address is out of bounds, returns -1;
	int getCube(int x, int y, int z);

	///dim is 0, 1, or 2, for x, y, z.  val is an integer.  if val is >= 0 and < the dimension of
	//the selected dimension, return true.  else false.
	bool inRange(int dim, int val);

	///returns the cubes which this triangle intersects
	set_t getCubes(double * t);

	///returns the cubes which this sphere intersects
	set_t getCubes(double * pt, double rad);

	///build the hash
	void build();

	/////////////////////////////////////////////////////////////////////////////
	///this gets the list of latticeObjs u want for a triangle insertion (Called by build)
	set_t getBoundingBoxOfObjs(double * t);
	/////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////////////////////
	///this finds the distance to the closest triangle(s) (this is a query)
	///tindex[0] is the triangle found.
	set_t distTo(double * pt, double * outputDist);
	/////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////
	///this determines if the given point is inside the surface
	///dist[0] is HUGE_VAL if outside the cube, otherwise it gives actual range
	///*****This method, which uses triangle normals, is deprecated.*****
	///bool isInside(double * pt, double * dist);
	///this determines if the given point is inside the surface
	///this is accomplished by casting a long ray (Segment) from pt, and counting the number
	/// of logical intersections (identical intersection points are counted once).  Even
	//  numbers of intersections mean that pt is outside, odd means that pt is inside.
	///USED TO BE CALLED THIS: bool isInsideRayTest(double * pt, double * distanceOutput);
	bool isInsideRayTest(double * pt, double * dist);
	bool isInsideRayTestMajority(double * pt, int numTests, bool * unanimous);
	/////////////////////////////////////////////////////////////////////////////
	//double isInsideWinding(double * pt);
	int isInsideClosestTriangle(double * pt);
	///this is the wrapper for insideOutside testing.
	bool isInside(double * pt);
	bool isInside(double ptx, double pty, double ptz);
	///this is a cache that stores old values for our repeated testing.
	int isInsideCached(double * pt);
	void ptCacheDeposit(double * pt, bool ioVal);
	void setCacheSize(int s);
	
	///This func returns -1 if outside, +1 if inside, 0 if failure.
	int isInsideLocalCubes(double * pt);
	
	///This returns 0 if inside the bounds, otherwise, -1 if not inside out, 1 if inside out.
	//int isinsideBoundsCheck(double * pt);

	/////////////////////////////////////////////////////////////////////////////
	///using flood fill, marks all cubes not containing geometry (Triangles) by setting their highlights to true;
	void floodFillMarker( set_t latticeIndices, int highlightCode );
	void floodFillStart(double * coord, int highlightCode);
	void floodFillStart(double x, double y, double z, int highlightCode);
	///returns true if inside a highlighted cube(s).  false otherwise.
	int insideHighlightedCube(double * coord);
	///eliminates all cube highlights (sets to false)
	void clearAllHighlights();
	/////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////
	///find Centroid of intersection points between input segment and the surface
	double * getAverageIntersectionPt(double * a, double * b);
	double * getAverageIntersectionPtNoTest(double * a, double * b);
	///counts the number of times a segment intersects with triangles from the surface
	int countIntersectionPts(double * a, double * b, double * dist);
	/////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////////////////////
	///Benchmarking function counts the number of flagged points, unflagged points.	
	void surveyFlagStatus();
	///Benchmarking function counts how many calls were made from external query.
	void surveyTestStatus();
	/////////////////////////////////////////////////////////////////////////////
	
	
	/////////////////////////////////////////////////////////////////////////////
	//This function tells you if a given triangle intersects with a triangle in the surface
	bool triangleIntersection( double * ta, double * tb, double * tc );	
	
	

#define TRIHASHCACHESIZE 10;

	SurfaceObject * surf;	//the Surface Data
	double res; 			//the size of each side of each cube.
	int *** grid;    		//a set of latticeObj indices (nesting all dimensions)
	LatticeObj ** objs;		//the actual objects are here.
	
	int xdim, ydim, zdim; //the dimension of each side of the grid.
	double xmin, xmax;	//extents for the cube.
	double ymin, ymax;
	double zmin, zmax;

	char * name;	///file that this was generated from (Debugging purposes)
	
	double * ioCacheVecs;	///cache of coordinates of tested values.
	bool * ioCacheBools;	///cache of inside/outside test results;
	int oldCachePtr;		///this points to the oldest entry.
	int newCachePtr;		///this points to the newest entry.
	int inCacheSize;		///the size of the cache. - the cache wraps around.
	
	bool isInsideOut; 		///This is true if the surface is to be evaluated in reverse.
						///(isInside() returns the opposite of what it should.)
						///This defaults to false in the constructor; you must set it to true;
						///added to implement negative spaces

	bool isEmpty;			///True if empty, false otherwise.
	int num_cache_hits;
	int num_local_cubes_tests;
	int num_long_ray_tests;
	int num_highlit_cube_tests;
};

#endif




