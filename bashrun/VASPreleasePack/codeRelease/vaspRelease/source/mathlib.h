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
 * File: mathlib.cpp
 *       interface for basic math library
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


#if !defined(_MATHLIB_H_)
#define _MATHLIB_H_


#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

double * normalizeVector(double * vector);
double * normalizeVector(double v1, double v2, double v3);
double vectorSize(double x, double y, double z);
double vectorSize(double * A, double * B);
double * crossProd(double * first, double * second);
double * getCentroid(double * vertices, int size);



	///wrapper added by Brian Chen, 5/15/2007
	///	note that order of corner enumeration is important!!!!
	///	It is assumed that the four corners make a circuit around the quad 
	///	i.e. q1, q2, q3, q4 draws the quad, not a "bow-tie"
	bool triQuadOverlapTest( double tcoords[9], double q0[3], double q1[3], double q2[3], double q3[3] );

	///wrapper added by Brian Chen, 8/03/2007
	///	note that order of corner enumeration is important!!!!
	///	It is assumed that the four corners make a circuit around the quad 
	///	i.e. q1, q2, q3, q4 draws the quad, not a "bow-tie"
	bool segQuadOverlapTest( double * s1, double * s2, double q0[3], double q1[3], double q2[3], double q3[3] );


	///OLD VERSION DEPRECATED - SEE CPP FILE
	//	double * intersect_Seg_Triangle( double * r0, double * r1, double * t0, double * t1, double * t2);
	//	double * intersect_Seg_Triangle_Slow( double * r0, double * r1, double * t0, double * t1, double * t2);
	///This was coded with assistance from http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/
	///Segment-Triangle Intersection
	double * intersect_Seg_Triangle( double * r0, double * r1, double * t0, double * t1, double * t2);
	
	///segment-Sphere Intersection
	///Coded with assistance from http://local.wasp.uwa.edu.au/~pbourke/geometry/sphereline/
	double * intersect_Seg_Sphere( double * r0, double * r1, double * c, double r);
	
	///finds the closest intersecting point from the segments from this triangle
	double * intersect_Triangle_Segments(double * t1, double * t2, double * t3, double * s0, double * s1);

	bool pointSegmentOverlapTest(double * pt, double * s1, double * s2);

	bool isColinear(double * pt1, double * pt2, double * pt3);
	bool isColinear(double * coords, int a, int b, int c);

	///turns a colinear triangle into the longest segment 
	///(returns 0 if pt2,pt3 is the longest, 1, if pt1,pt3 is the longest, 2 if pt1,pt2 is the longest)
	///assumes that the triangle is actually colinear; otherwise returns longest edge.
//	int extractSegment(double * pt1, double * pt2, double * pt3);

	///tests if a segment AB intersects an axis aligned quad, defined by ranges x1, x2, y1, y2, z1, z2.
	///if intersect is non-NULL, intersect[0], interect[1], intersect[2] is filled with the intersection pt.
	/// returns 1 if intersects, 0 if not intersect.
	int segAAquadCollision(double * A, double * B,
						double x1, double x2, double y1, double y2, double z1, double z2, 
						double * intersect );

	///determines if a triangle intersects with an axis aligned cube.
	bool triCubeIntersect(double * tcoords, double xn, double xp, double yn, double yp, double zn, double zp);






/*----------------Macros----------------*/
///Cross Product
#ifndef CROSS
#define CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#endif

///Dot Product 
#ifndef DOT
#define DOT(v1,v2) ((v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]))
#endif



#endif


/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/




