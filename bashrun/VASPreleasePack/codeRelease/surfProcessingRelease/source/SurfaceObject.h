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
 * File: SurfaceObject.cpp
 *       interface for Object for representing triangle meshes
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

#ifndef SURFACEOBJECT_H
#define SURFACEOBJECT_H

#include "StdAfx.h"

class SurfaceObject  
{
public:
	int numPoints;
	int numTriangles;
	
	///mandatory memory allocation
	double * surfacePoints;
	double * surfaceNormals;		///these averaged normals are used for rendering
	int * triangles;
	double * triangleNormals;	///these normals are used for flat shading and interior checking
	double * centroid;

	///optional memory allocation
	set_t highlights;
	set_t edges; 				///edge highlights
	double * colors;			///optional colors, initialized only through separate request.  One vector for each pt.


	///////////////////////////////////////////////////////
	///Constructors / Destructors
	SurfaceObject();
	SurfaceObject(int numPts, double * ptsAndNorms, int numTris, int * inputTris, bool elimIntCavities);
	SurfaceObject(int numPts, double * pts, double * norms, int numTris, int * tris, double * triangleNorms);
	SurfaceObject( set_t surfSet );
	virtual ~SurfaceObject();
	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////
	///returns the coords of the requested triangle
	///returns null if out of range.
	///notice that the triangle must be formatted in this way in the array (3 then 3 then 3)
	double * getTriangle(int num); 
	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////
	///Computes the centroid based on the point positions
	double * getCentroid(); 
	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////
	///Adds the geometry of another SurfaceObject to this SurfaceObject.
	///ASSUMES THAT THIS SURFACEOBJECT DOES NOT COLLIDE WITH THE OTHER SURFACE OBJECT.
	///DOES NOT CHECK FOR COLLISION, SO MUST BE DONE BEFORE THIS FUNCTION IS CALLED.
	void addObject(SurfaceObject * obj);
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///Copies this surfaceObject into another surfaceObject
	///If indices is non-NULL, returns a subset of the surface
	///using only the points that are the indices of this surface.
	SurfaceObject * copy(set_t indices);
	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////
	///identifies edges which are not surrounded by triangles - 
	///locations where the surface is non manifold.
	///This turns the edges variable non-null.
	///
	///THIS FUNCTION IS CURRENTLY BROKEN
	//
	void identifyNonManifoldEdges(void);
	///helper function tells you how int elements are identical between two sets.
	int countNumCommon(set_t set1, set_t set2);
	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////
	///This function flips the normals of the surfaceObject backwards 
	///so that we can treat it as a negative volume
	void flipNormals();
	///////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////
	///add colors
	void addColors(double * c);
	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////
	///print out details
	void toString();
	void printSummary();
	///////////////////////////////////////////////////////

};

#endif 



