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
 * File: LatticeObj.cpp
 *       Local information for lattice hashing (interface)
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

/*
##############################################################################################
# ####  ##   #   #   ##  ##   ###  ##   ##    #   ##  ##   ###################################
# ### ## ## ### ### ## ## # #### ## # ## ### ## ### ## ## ####################################
# ### ## ## ### ### ## ####  ### ## #   #### ##  ## ##### ####################################
# ###    ## ### ### ## #### #### ## # ## ### ## ### ##### ####################################
# ### ## ## ### ### ## ## # #### ## # ## # # ## ### ## ## ####################################
#   # ## ## ### ##   ##  ##   ###  ##   ##   ##   ##  ### ####################################
##############################################################################################
*/

#ifndef LATTICE_OBJECT_H
#define LATTICE_OBJECT_H

#include "StdAfx.h"


#define HIGHLIGHT_INTERIOR 100
#define HIGHLIGHT_EXTERIOR 200
#define NOT_HIGHLIGHTED    000

class SurfaceObject;

class LatticeObj
{
	public:

	///constructor for latticeObj
	///  cors fills in for corners
	///  ind is the index of thsi latticeObj
	///  a is the array of adjancies (always 26, subbing -1 for non-adjacencies.
	///  sobj is the surfaceObj
//	LatticeObj(  );
	LatticeObj( double * cors, int ind, SurfaceObject * sobj );
	virtual ~LatticeObj();

	///set the adjacency data;
	///  a is the array of adjancies (always 27), subbing -1 for non-adjacencies.
	///  one is the reflexibe adjacency
	void setAdjacency(int * a);

	///set the triangle data;
	void addTriangle(int t);
	
	///determine if a triangle intersects the volume of the cube.
	bool triangleIntersect( double * tcoords );
	

	///determine if a segment intersects the volume of the cube.
	bool segmentIntersect( double * a, double * b );
	///This old apparatus was developed to test a new implementation of segmentIntersect
	///and is now deprecated.
//	bool segmentIntersect( double * a, double * b );
//	bool segmentIntersectOld( double * a, double * b );
	
	///Finds the smallest distance to all triangles in the cube to the point
	///returns HUGE_VAL if cube is empty.
	///tIndex is the triangle index of the closest triangle
	double getClosestDist( double * p, int &tIndex );  
	
	///Finds the intersection Pts with all triangles that intersect segment a-b
	//output array is packed by first having the size, then 3xVectors
	//Input arrays are 3x vectors
	double * getIntersectPts( double * a, double * b );

	int index;			///the global ID of this cube
	double * corners;	///the corners of this cube. (8*3 doubles, for queries and insertion)
	set_t adj;		///the indices of adjacent cubes.	(this is for expanding queries)
	set_t indices;		///the triangle indices we will store here (this is to answer queries)
//	double maxX, minX, maxY, minY, maxZ, minZ;	///extents of the cube.
	SurfaceObject * surf; 	///utility reference to the SurfaceObject for testing.
	int highlight;


};


///determines if a triangle intersects with an axis aligned cube.
bool triCubeIntersect(double * tcoords, double xn, double xp, double yn, double yp, double zn, double zp);


#endif

