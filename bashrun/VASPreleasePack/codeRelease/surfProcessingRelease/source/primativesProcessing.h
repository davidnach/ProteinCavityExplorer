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
 * File: primativesProcessing.h
 *       interface for surface generation for abstract spheres.
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

#ifndef PRIMATIVE_PROCESSING
#define PRIMATIVE_PROCESSING

#define PRIM_TYPE_SPHERE 100
#define PRIM_TYPE_CUBE 101

class PrimSphere;
class PrimCube;
class LatticeHash;
class TriangleLatticeHash;

////////////////////////////////////////////////////////
// Primative class is a general class for all CSG primatives
// Handles all intersection cases.
//
class Primative
{
public:
	Primative();
	virtual ~Primative();

	int primType;
	
	virtual bool isInside(double * pt);
	///it is understood that pt1 is always inside
	virtual double * intersectionPt(double * pt1, double * pt2);
	virtual double * getBounds();
	virtual Primative * copy();
	virtual SurfaceObject * genSurface();
	virtual double * normal(double * pt);
	virtual double distance(double * pt);

	bool intersects( Primative * p );
	bool sphereIntersectsSphere( PrimSphere * s1, PrimSphere * s2);	
	bool sphereIntersectsCube  ( PrimSphere * s1, PrimCube * c2);	
	bool cubeIntersectsCube    ( PrimCube   * c1, PrimCube * c2);	
};




////////////////////////////////////////////////////////
// Primative Sphere class
//
class PrimSphere : Primative
{
public:
	PrimSphere(double * o, double r);
	virtual ~PrimSphere();
	
	double * origin;	///double[3], x, y, z
	double radius;		///double radius of sphere
	
	double * getBounds();
	bool isInside(double * pt);
	double * intersectionPt(double * pt1, double * pt2);
	PrimSphere * copy();
	SurfaceObject * genSurface();
	double * normal(double * pt);
	double distance(double * pt);

};



////////////////////////////////////////////////////////
// Primative cube class
//
class PrimCube : Primative
{
public:
	PrimCube(double * c, double s);
	virtual ~PrimCube();
	
	double * center;	///double[3], x, y, z - center of cube
	double size;		///one half the sidelength of the cube.
	
	double * getBounds();
	bool isInside(double * pt);
	double * intersectionPt(double * pt1, double * pt2);
	PrimCube * copy();
	SurfaceObject * genSurface();
	double * normal(double * pt);
	double distance(double * pt);

};



////////////////////////////////////////////////////////
// Container class for multiple primatives
//
class primativesSet
{
public:
	primativesSet();
	primativesSet(char * filename);
	virtual ~primativesSet();

	///////////////////////////////
	///Variables
	///////////////////////////////
	set_t objs;
	bool isInsideOut;			//true if this is inside out, false otherwise.
	LatticeHash * hash;			//NULL unless constructHashTable is called.
	
	///////////////////////////////
	///Methods 
	///////////////////////////////
	int size();				///returns number of objs in objs
	Primative * getPrim(int ind);		///returns a pointer to the obj desired; NULL if not in range.
	void addObj(Primative * obj);

	///tests if a point is inside the primatives:
	bool isInside(double * pt);			//dispatch.
	bool isInsideLinear(double * pt);		//default; checks all primatives for every point.
	bool isInsideHashed(double * pt);		//uses a hash table to check only nearby primatives.
	
	///it is understood that pt1 is always inside
	double * intersectionPt(double * pt1, double * pt2);		//dispatch
	double * intersectionPtLinear(double * pt1, double * pt2);	//default; checks all primatives for every point.
	double * intersectionPtHashed(double * pt1, double * pt2);	//uses a hash table to check only nearby primatives




	double * getPrimSetBounds();
	double * normal(double * pt);	///find the closest primative, then compute the normal
	
	void constructHashTable();	//sets up the LatticeHash

	void toString();

//	void outputPrimsFile(char * fileName);
	
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//these functions turn probesets into SurfaceObjects, and know how to change a simple primative
///into a simple object, rather than using heavy marching cubes.

///wrapper parser calls probeGen after generating a set of spheres from a pdb file
SurfaceObject * ProbeGenWrapper(char * ligFile, double Sphere_Radius, double resolution);
///wrapper parser calls probeGen after generating a set of primatives from a PRIMS file
SurfaceObject * ProbeGenWrapper(char * primsFile, double resolution);


///this function identifies the separate stuff, and generates those with simple geometry.
SurfaceObject * convertPrimativesToSurface( primativesSet * primSet, double resolution );

//this function identifies the hard stuff and generates that with marching cubes (probeGen)
SurfaceObject * convertPrimsMarchingCubes( primativesSet * primSet, double resolution );

//Generates a surface from a PrimitivesSet using marching cubes.
SurfaceObject * ProbeGen( primativesSet * primSet, double resolution );


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
///Boolean marching cubes using a probeset.  Leverages the accuracy of the 
///probeset's analytic computations.

///custom boolean Geometry processing for heterogeneous data types - one hash one primSet
double * getBooleanPrimGeometry(bool p1a, bool p1b, bool p2a, bool p2b, 
						TriangleLatticeHash * hasha, primativesSet * primSet, 
						double * p1, double * p2, int op);

///helper function adapted for primatives processing.
SurfaceObject * booleanProbeProcessing(SurfaceObject * surf1, bool neg1, primativesSet * primSet, bool neg2, int op, double res);







#endif




