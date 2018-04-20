/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 *   Copyright (c) 2014 Brian Y. Chen
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
 * This file incorporates work covered by the following copyright
 * and permission notices:
 *
 *   Copyright (c) 2001 softsurfer.com
 *   Copyright (c) 2002, 2003 Philippe Guigue and Olivier Devillers
 *
 * Permission to use, copy, modify, and/or distribute this software
 * for any purpose with or without fee is hereby granted, provided
 * that the above copyright notice and this permission notice appear
 * in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
 * AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFTS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * File: triIntersect.cpp
 *       geometric tests for triangles
 *
 * Written by 
 *       softSurfer.com
 *       Philippe Guigue, Olivier Devillers
 *
 ***************************************************************************/


#include "triIntersect.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This copyright applies to the code below, until the identical statement of this copyright.  It has been modified by Brian Chen
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
References
David Eberly, 3D Game Engine Design, Section 2.6 "Distance Methods" (2001)
Seth Teller, line_line_closest_points3d() (2000) cited in the Graphics Algorithms FAQ Subject 5.18 (2001)
*/
//===================================================================
// dist3D_Segment_to_Segment():
//    Input:  two 3D line segments S1 and S2
//    Return: the shortest distance between S1 and S2
//    Secondary Output: closest - this is the point on s2 that is closest to s1
//                      if closest is NULL, nothing is returned for closest, without segfault.
//
//  Modified for usage in VASP, by Brian Chen, 2010.
///////////////////////////////////////////////////////////////////////////////////
double dist3D_Segment_to_Segment( double * s1p0, double * s1p1, double * s2p0, double * s2p1, double * closest)
{
	//Vector   u = S1.P1 - S1.P0;
	double * u = new double[3];
	u[0] = s1p1[0] - s1p0[0];
	u[1] = s1p1[1] - s1p0[1];
	u[2] = s1p1[2] - s1p0[2];
	
	//Vector   v = S2.P1 - S2.P0;
	double * v = new double[3];
	v[0] = s2p1[0] - s2p0[0];
	v[1] = s2p1[1] - s2p0[1];
	v[2] = s2p1[2] - s2p0[2];
	
	//Vector   w = S1.P0 - S2.P0;
	double * w = new double[3];
	w[0] = s1p0[0] - s2p0[0];
	w[1] = s1p0[1] - s2p0[1];
	w[2] = s1p0[2] - s2p0[2];
	
	double    a = DOT(u,u);        // always >= 0
	double    b = DOT(u,v);
	double    c = DOT(v,v);        // always >= 0
	double    d = DOT(u,w);
	double    e = DOT(v,w);
	double    D = a*c - b*b;       // always >= 0
	double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
	double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < SMALL_NUM) {  // the lines are almost parallel
		sN = 0.0;        // force using point P0 on segment S1
		sD = 1.0;        // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	// get the closest points on the infinite lines
	else {
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else{ 
			if (sN > sD) {  // sc > 1 => the s=1 edge is visible
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}
	}
	
	if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0){
			sN = 0.0;
		}
		else{ 
			if (-d > a){
				sN = sD;
			}
			else{
				sN = -d;
				sD = a;
			}
		}
	}
	else{
		// tc > 1 => the t=1 edge is visible
		if (tN > tD) {
			tN = tD;
			// recompute sc for this edge
			if ((-d + b) < 0.0){
				sN = 0;
			}
			else{
				if ( (-d + b) > a ){
					sN = sD;
				}
				else {
					sN = (-d + b);
					sD = a;
				}
			}
		}
	}

	// finally do the division to get sc and tc
	if( fabs(sN) < SMALL_NUM ){ sc = 0.0; }
	else{ sc = (sN/sD); }

	if( fabs(tN) < SMALL_NUM ){ tc = 0.0; }
	else{ tc = (tN/tD); }

	// get the difference of the two closest points
	//Vector   dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
	double * dP = new double[3];
	dP[0] = w[0] + sc*u[0] - tc*v[0];
	dP[1] = w[1] + sc*u[1] - tc*v[1];
	dP[2] = w[2] + sc*u[2] - tc*v[2];

	// Now we assign the closest point
	if(closest != NULL){
		closest[0] = tc*v[0];
		closest[1] = tc*v[1];
		closest[2] = tc*v[2];
	}

	// return the closest distance
	double result = vectorSize( dP[0], dP[1], dP[2] );

	delete[](u);
	delete[](v);
	delete[](w);
	delete[](dP);

	return result;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This copyright applies to the code above, up until the identical statement of this copyright.  It has been modified by Brian Chen
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//The remainder of the code in this file is cited by this work below.  It was written by P. Guigue and O. Devillers 
/*
*		  		
*  Triangle-Triangle Overlap Test Routines				
*  July, 2002                                                          
*  Updated December 2003                                                
*                                                                       
*  This file contains C implementation of algorithms for                
*  performing two and three-dimensional triangle-triangle intersection test 
*  The algorithms and underlying theory are described in                    
*                                                                           
* "Fast and Robust Triangle-Triangle Overlap Test 
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*                                                 
*  Journal of Graphics Tools, 8(1), 2003                                    
*                                                                           
*  Several geometric predicates are defined.  Their parameters are all      
*  points.  Each point is an array of two or three double precision         
*  floating point numbers. The geometric predicates implemented in          
*  this file are:                                                            
*                                                                           
*    int tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)                         
*    int tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)                         
*                                                                           
*    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
*                                     coplanar,source,target)               
*                                                                           
*       is a version that computes the segment of intersection when            
*       the triangles overlap (and are not coplanar)                        
*                                                                           
*    each function returns 1 if the triangles (including their              
*    boundary) intersect, otherwise 0                                       
*                                                                           
*                                                                           
*  Other information are available from the Web page                        
*  http://www.acm.org/jgt/papers/GuigueDevillers03/                         
*                                                                           
*/



///////////////////////////////////////////////////////////////////////////////////
// Three-dimensional Triangle-Triangle Overlap Test
///////////////////////////////////////////////////////////////////////////////////
int tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3], double p2[3], double q2[3], double r2[3])
{
  double dp1, dq1, dr1, dp2, dq2, dr2;
  double v1[3], v2[3];
  double N1[3], N2[3]; 
  
  /* Compute distance signs  of p1, q1 and r1 to the plane of
     triangle(p2,q2,r2) */

  SUB(v1,p2,r2)
  SUB(v2,q2,r2)
  CROSS(N2,v1,v2)

  SUB(v1,p1,r2)
  dp1 = DOT(v1,N2);
  SUB(v1,q1,r2)
  dq1 = DOT(v1,N2);
  SUB(v1,r1,r2)
  dr1 = DOT(v1,N2);
  
  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0; 

  /* Compute distance signs  of p2, q2 and r2 to the plane of
     triangle(p1,q1,r1) */
  
  SUB(v1,q1,p1)
  SUB(v2,r1,p1)
  CROSS(N1,v1,v2)

  SUB(v1,p2,r1)
  dp2 = DOT(v1,N1);
  SUB(v1,q2,r1)
  dq2 = DOT(v1,N1);
  SUB(v1,r2,r1)
  dr2 = DOT(v1,N1);
  
  if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

  /* Permutation in a canonical form of T1's vertices */

  if (dp1 > 0.0f) {
    if (dq1 > 0.0f) TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)	
    else TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } else if (dp1 < 0.0f) {
    if (dq1 < 0.0f) TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } else {
    if (dq1 < 0.0f) {
      if (dr1 >= 0.0f) TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) {
      if (dr1 > 0.0f) TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  {
      if (dr1 > 0.0f) TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
    }
  }
}


int coplanar_tri_tri3d(double p1[3], double q1[3], double r1[3], double p2[3], double q2[3], double r2[3], double normal_1[3], double normal_2[3]){
  
  double P1[2],Q1[2],R1[2];
  double P2[2],Q2[2],R2[2];

  double n_x, n_y, n_z;

  n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
  n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
  n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


  /* Projection of the triangles in 3D onto 2D such that the area of
     the projection is maximized. */

  if (( n_x > n_z ) && ( n_x >= n_y )) {
    // Project onto plane YZ

      P1[0] = q1[2]; P1[1] = q1[1];
      Q1[0] = p1[2]; Q1[1] = p1[1];
      R1[0] = r1[2]; R1[1] = r1[1]; 
    
      P2[0] = q2[2]; P2[1] = q2[1];
      Q2[0] = p2[2]; Q2[1] = p2[1];
      R2[0] = r2[2]; R2[1] = r2[1]; 

  } else if (( n_y > n_z ) && ( n_y >= n_x )) {
    // Project onto plane XZ

    P1[0] = q1[0]; P1[1] = q1[2];
    Q1[0] = p1[0]; Q1[1] = p1[2];
    R1[0] = r1[0]; R1[1] = r1[2]; 
 
    P2[0] = q2[0]; P2[1] = q2[2];
    Q2[0] = p2[0]; Q2[1] = p2[2];
    R2[0] = r2[0]; R2[1] = r2[2]; 
    
  } else {
    // Project onto plane XY

    P1[0] = p1[0]; P1[1] = p1[1]; 
    Q1[0] = q1[0]; Q1[1] = q1[1]; 
    R1[0] = r1[0]; R1[1] = r1[1]; 
    
    P2[0] = p2[0]; P2[1] = p2[1]; 
    Q2[0] = q2[0]; Q2[1] = q2[1]; 
    R2[0] = r2[0]; R2[1] = r2[1]; 
  }

  return tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);
    
}


///////////////////////////////////////////////////////////////////////////////////
// Three-dimensional Triangle-Triangle Intersection              
//
// The following version computes the segment of intersection of the
// two triangles if it exists. 
// coplanar returns whether the triangles are coplanar
// source and target are the endpoints of the line segment of intersection 
///////////////////////////////////////////////////////////////////////////////////
int tri_tri_intersection_test_3d(double p1[3], double q1[3], double r1[3],  double p2[3], double q2[3], double r2[3], int * coplanar, double source[3], double target[3] )
				 
{
  double dp1, dq1, dr1, dp2, dq2, dr2;
  double v1[3], v2[3], v[3];
  double N1[3], N2[3], N[3];
  double alpha;

  // Compute distance signs  of p1, q1 and r1 
  // to the plane of triangle(p2,q2,r2)

  SUB(v1,p2,r2)
  SUB(v2,q2,r2)
  CROSS(N2,v1,v2)

  SUB(v1,p1,r2)
  dp1 = DOT(v1,N2);
  SUB(v1,q1,r2)
  dq1 = DOT(v1,N2);
  SUB(v1,r1,r2)
  dr1 = DOT(v1,N2);
  
  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0; 

  // Compute distance signs  of p2, q2 and r2 
  // to the plane of triangle(p1,q1,r1)
  
  SUB(v1,q1,p1)
  SUB(v2,r1,p1)
  CROSS(N1,v1,v2)

  SUB(v1,p2,r1)
  dp2 = DOT(v1,N1);
  SUB(v1,q2,r1)
  dq2 = DOT(v1,N1);
  SUB(v1,r2,r1)
  dr2 = DOT(v1,N1);
  
  if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

  // Permutation in a canonical form of T1's vertices

  if (dp1 > 0.0f) {
    if (dq1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
	
    else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } else if (dp1 < 0.0f) {
    if (dq1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } else {
    if (dq1 < 0.0f) {
      if (dr1 >= 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) {
      if (dr1 > 0.0f) TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  {
      if (dr1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else {
	// triangles are co-planar

	*coplanar = 1;
	return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////
// Two dimensional Triangle-Triangle Overlap Test    
///////////////////////////////////////////////////////////////////////////////////
int ccw_tri_tri_intersection_2d(double p1[2], double q1[2], double r1[2], 
				double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
      else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
}


int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
			    double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p1,q1,r1) < 0.0f )
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/





