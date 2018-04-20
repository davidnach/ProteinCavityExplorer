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
 *       basic math library
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

////////////////////
//Returns normalized vector in input direction
double * normalizeVector(double * vector)
{
	double * result = new double[3];

	////calculate the magnitude
		double mag = sqrt(	(vector[0]*vector[0]) + 
							(vector[1]*vector[1]) + 
							(vector[2]*vector[2]) );

	////set the normalized vector
		result[0] = vector[0]/mag;
		result[1] = vector[1]/mag;
		result[2] = vector[2]/mag;

	return result;

}

double * normalizeVector(double v0, double v1, double v2)
{
	double * result = new double[3];
	double mag = sqrt( (v0*v0) + (v1*v1) + (v2*v2) );
	result[0] = v0/mag;
	result[1] = v1/mag;
	result[2] = v2/mag;

	return result;	
}

///finds the matrix magnitude
double vectorSize(double x, double y, double z)
{
	return sqrt( (x*x) + (y*y) + (z*z) );
}

///finds the matrix magnitude
double vectorSize(double * A, double * B)
{
	double x = B[0]-A[0];
	double y = B[1]-A[1];
	double z = B[2]-A[2];

	return sqrt( (x*x) + (y*y) + (z*z) );
}



double * crossProd(double * first, double * second)
{
	///Fix size of array to 3 doubles
	double * result = new double[3];	

	/////calculate cross product
		result[0] = (first[1] * second[2]) - (first[2] * second[1]);
		result[1] = -( (first[0] * second[2]) - (first[2] * second[0]) );
		result[2] = (first[0] * second[1]) - (first[1] * second[0]);

	return result;

}

double * getCentroid(double * vertices, int size)
{
	double * result = new double[3];
	double tempSum = 0.0;
	int i = 0;
	
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+0];}
	result[0] = tempSum/size;

	tempSum = 0.0;
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+1];}
	result[1] = tempSum/size;

	tempSum = 0.0;
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+2];}
	result[2] = tempSum/size;

	return result;
}






///////////////////////////////////////////////////////////////////////////////////
///wrapper added by Brian Chen, 8/03/2007
///	note that order of corner enumeration is important!!!!
///	It is assumed that the four corners make a circuit around the quad 
///	i.e. q1, q2, q3, q4 draws the quad, not a "bow-tie"
///////////////////////////////////////////////////////////////////////////////////
bool segQuadOverlapTest( double * s1, double * s2, double q0[3], double q1[3], double q2[3], double q3[3] )
{
	bool result = false;
	
	///test overlaps.  make sure triangles are right.
	double * info1 = intersect_Seg_Triangle( s1, s2, q0, q1, q2);
	double * info2 = intersect_Seg_Triangle( s1, s2, q2, q3, q0);

	///test if either overlaps
	if(info1!=NULL){ result = true; delete[](info1); }
	if(info2!=NULL){ result = true; delete[](info2); }
	return result;
}


///////////////////////////////////////////////////////////////////////////////////
///finds the closest intersecting point from the segments from this triangle
//by Brian Chen
///////////////////////////////////////////////////////////////////////////////////
double * intersect_Triangle_Segments(double * t1, double * t2, double * t3, double * s0, double * s1)
{
	double * s1p1 = new double[3];	s1p1[0] = t1[0];	s1p1[1] = t1[1];	s1p1[2] = t1[2];
	double * s1p2 = new double[3];	s1p2[0] = t2[0];	s1p2[1] = t2[1];	s1p2[2] = t2[2];
	double * s2p1 = new double[3];	s2p1[0] = t2[0];	s2p1[1] = t2[1];	s2p1[2] = t2[2];
	double * s2p2 = new double[3];	s2p2[0] = t3[0];	s2p2[1] = t3[1];	s2p2[2] = t3[2];
	double * s3p1 = new double[3];	s3p1[0] = t3[0];	s3p1[1] = t3[1];	s3p1[2] = t3[2];
	double * s3p2 = new double[3];	s3p2[0] = t1[0];	s3p2[1] = t1[1];	s3p2[2] = t1[2];

/*
	printf("s1p1: %f %f %f\n", s1p1[0], s1p1[1], s1p1[2]);
	printf("s1p2: %f %f %f\n", s1p2[0], s1p2[1], s1p2[2]);
	printf("s2p1: %f %f %f\n", s2p1[0], s2p1[1], s2p1[2]);
	printf("s2p2: %f %f %f\n", s2p2[0], s2p2[1], s2p2[2]);
	printf("s3p1: %f %f %f\n", s3p1[0], s3p1[1], s3p1[2]);
	printf("s3p2: %f %f %f\n", s3p2[0], s3p2[1], s3p2[2]);
*/

	double * closestPt = new double[3];
	double * closePt = new double[3];
	double dist;
	double minDist = HUGE_VAL;
	
	if(vectorSize(s1p1, s1p2) > 0){    ///check for positive distance, to avoid degeneracy.
		dist = dist3D_Segment_to_Segment( s1p1, s1p2, s0, s1, closePt);
		if(dist < minDist){ closestPt[0] = closePt[0]; closestPt[1] = closePt[1]; closestPt[2] = closePt[2]; minDist = dist;}
	}
//	printf("TEST DISTANCE1: %f\n", minDist);
	if(vectorSize(s2p1, s2p2) > 0){    ///check for positive distance, to avoid degeneracy.
		dist = dist3D_Segment_to_Segment( s2p1, s2p2, s0, s1, closePt);
		if(dist < minDist){ closestPt[0] = closePt[0]; closestPt[1] = closePt[1]; closestPt[2] = closePt[2]; minDist = dist;}
	}
//	printf("TEST DISTANCE2: %f\n", minDist);
	if(vectorSize(s3p1, s3p2) > 0){    ///check for positive distance, to avoid degeneracy.
		dist = dist3D_Segment_to_Segment( s3p1, s3p2, s0, s1, closePt);
		if(dist < minDist){ closestPt[0] = closePt[0]; closestPt[1] = closePt[1]; closestPt[2] = closePt[2]; minDist = dist;}
	}
//	printf("TEST DISTANCE3: %f\n", minDist);

	
	delete[](s1p1); delete[](s1p2);
	delete[](s2p1); delete[](s2p2);
	delete[](s3p1); delete[](s3p2);
	delete[](closePt);

	if( minDist < SMALL_NUM ){
		return closestPt;
	}
	delete[](closestPt);
	return NULL;

}


bool pointSegmentOverlapTest(double * pt, double * s1, double * s2)
{
	double * t1 = new double[3];
	t1[0] = pt[0]-s1[0];	t1[1] = pt[1]-s1[1];	t1[2] = pt[2]-s1[2];
	double * t2 = new double[3];
	t2[0] = s2[0]-s1[0];	t2[1] = s2[1]-s1[1];	t2[2] = s2[2]-s1[2];

	bool result = false;
	double dprod = DOT(t1, t2);

	delete[](t1);	delete[](t2);
	if(dprod > 0 && dprod < 1){ result = true; }

	return result;
}



///////////////////////////////////////////////////////////////////////////////////
// if these three points are fairly collinear, then return true; otherwise false.
///////////////////////////////////////////////////////////////////////////////////
bool isColinear(double * pt1, double * pt2, double * pt3)
{
	bool result = false;
	
	double d1 = vectorSize(pt1, pt2);
	double d2 = vectorSize(pt2, pt3);
	double d3 = vectorSize(pt3, pt1);

	int longestEdge = -1;
	double longest = 0;
	if(d1 > longest){ longestEdge = 1; longest = d1; }
	if(d2 > longest){ longestEdge = 2; longest = d2; }
	if(d3 > longest){ longestEdge = 3; longest = d3; }

	///if this is a point triangle, return true;
	if(longest < .001){
		result = true;  
		return result;
	}

	double otherSides = 0;
	if(longestEdge == 1){ otherSides = d2+d3; }
	if(longestEdge == 2){ otherSides = d3+d1; }
	if(longestEdge == 3){ otherSides = d1+d2; }

	if(otherSides - longest < .001){
		result = true;
	}

	return result;
}

bool isColinear(double * coords, int a, int b, int c)
{
	double * A = new double[3];
	A[0] = coords[3*a+0]; A[1] = coords[3*a+1]; A[2] = coords[3*a+2]; 
	double * B = new double[3];
	B[0] = coords[3*b+0]; B[1] = coords[3*b+1]; B[2] = coords[3*b+2]; 
	double * C = new double[3];
	C[0] = coords[3*c+0]; C[1] = coords[3*c+1]; C[2] = coords[3*c+2]; 

	bool result = isColinear(A, B, C);
	delete[](A);
	delete[](B);
	delete[](C);

	return result;
}




///////////////////////////////////////////////////////////////////////////////////
// Seg-Tri intersection.  Brian Chen, 2010
///////////////////////////////////////////////////////////////////////////////////
double * intersect_Seg_Triangle( double * r0, double * r1, double * t0, double * t1, double * t2)
{
	double * result = NULL;

	double ba0 = t1[0]-t0[0];	double ba1 = t1[1]-t0[1];	double ba2 = t1[2]-t0[2];
	double ca0 = t2[0]-t0[0];	double ca1 = t2[1]-t0[1];	double ca2 = t2[2]-t0[2];
	double norm0 = (ba1*ca2) - (ba2*ca1);
	double norm1 = -((ba0*ca2) - (ba2*ca0));
	double norm2 = (ba0*ca1) - (ba1*ca0);

	/////////////////////////////////////////////////////////////////////////
	///DEGENERACY CHECK FIRST:///////////////////////////////////////////////
	///deal with degenerate triangles: triangles that are essentially lines or points
	double check_area = vectorSize(norm0, norm1, norm2);
	double check_d1 = vectorSize(t0,t1);
	double check_d2 = vectorSize(t1,t2);
	double check_d3 = vectorSize(t2,t0);
	
	if(check_area < SMALL_NUM){
		//if this is a point-triangle, take any point of the triangle, see if it intersects.
		if(check_d1 < SMALL_NUM && check_d2 < SMALL_NUM && check_d3 < SMALL_NUM ){
			bool check_pt = pointSegmentOverlapTest(t0, r0, r1);
			if(check_pt){
				double * check_result = new double[3];
				check_result[0] = t0[0]; check_result[1] = t0[1]; check_result[2] = t0[2]; 
				return check_result;
			}
			else{ return result; }
		}
		//if this is a line-triangle, then
		else{
			result = intersect_Triangle_Segments(t0, t1, t2, r0, r1);
			return result;
		}
	}
	///DEGENERACY CHECK FIRST:///////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/*
	//	t1
	//   | \
	//	|  \
	//   u   \
	//	|    \
	//	t0-v-t2
	*/
	///////////////////////////////////////////////////////////////////////
	///first check if the segment even intersects the plane of the triangle.
	
	double mag;
	mag = sqrt( (ba0*ba0) + (ba1*ba1) + (ba2*ba2) );
	double u0 = ba0/mag; 	double u1 = ba1/mag; 	double u2 = ba2/mag; 
	mag = sqrt( (ca0*ca0) + (ca1*ca1) + (ca2*ca2) );
	double v0 = ca0/mag;	double v1 = ca1/mag;	double v2 = ca2/mag;	

	double uvCross0 = (u1*v2) - (u2*v1);
	double uvCross1 = -((u0*v2) - (u2*v0));
	double uvCross2 = (u0*v1) - (u1*v0);

	mag = sqrt( (uvCross0*uvCross0) + (uvCross1*uvCross1) + (uvCross2*uvCross2) );
	double tnorm0 = uvCross0/mag;
	double tnorm1 = uvCross1/mag;
	double tnorm2 = uvCross2/mag;

	mag = sqrt( ((r0[0]-t0[0])*(r0[0]-t0[0])) + ((r0[1]-t0[1])*(r0[1]-t0[1])) + ((r0[2]-t0[2])*(r0[2]-t0[2])) );
	if(mag==0){ double * result = new double[3]; result[0] = t0[0]; result[1] = t0[1]; result[0] = t0[2]; return result; }
	double r0t0 = (r0[0]-t0[0])/mag;	///vect0
	double r0t1 = (r0[1]-t0[1])/mag;
	double r0t2 = (r0[2]-t0[2])/mag;
	mag = sqrt( ((r1[0]-t0[0])*(r1[0]-t0[0])) + ((r1[1]-t0[1])*(r1[1]-t0[1])) + ((r1[2]-t0[2])*(r1[2]-t0[2])) );
	if(mag==0){ double * result = new double[3]; result[0] = t0[0]; result[1] = t0[1]; result[0] = t0[2]; return result; }
	double r1t0 = (r1[0]-t0[0])/mag;	///vect1
	double r1t1 = (r1[1]-t0[1])/mag;
	double r1t2 = (r1[2]-t0[2])/mag;

	///////Use the dot product to dtermine if the vector is with, against, or perpendicular to the normal
	double test0 = (tnorm0*r0t0) + (tnorm1*r0t1) + (tnorm2*r0t2);
	double test1 = (tnorm0*r1t0) + (tnorm1*r1t1) + (tnorm2*r1t2);

	///eliminate issues with small numbers
	if(fabs(test0) < SMALL_NUM){ test0 = 0.0; }
	if(fabs(test1) < SMALL_NUM){ test1 = 0.0; }

	///////if the dot product is 0, then it is perpendicular, and the point is on the plane.
	///////if the dot product is positive, then it is on the same side as the normal
	///////if the dot product is negative, then it is on the opposite side of the triangle from the normal.
	int triangleCase = -1;
	if( (test0<0)  && (test1<0)  ){ triangleCase = 0; }
	if( (test0<0)  && (test1==0) ){ triangleCase = 1; }
	if( (test0<0)  && (test1>0)  ){ triangleCase = 2; }
	if( (test0==0) && (test1<0)  ){ triangleCase = 3; }
	if( (test0==0) && (test1==0) ){ triangleCase = 4; }
	if( (test0==0) && (test1>0)  ){ triangleCase = 5; }
	if( (test0>0)  && (test1<0)  ){ triangleCase = 6; }
	if( (test0>0)  && (test1==0) ){ triangleCase = 7; }
	if( (test0>0)  && (test1>0)  ){ triangleCase = 8; }
		
	if(triangleCase == -1){ 
		printf("WTF!! triangleCase == -1 TEST0: [%f] TEST1: [%f]\n", test0, test1 ); 
		printf("ba0 %f  ba1 %f  ba2 %f  ca0 %f  ca1 %f  ca2 %f\n", ba0, ba1, ba2, ca0, ca1, ca2 );
		printf("TEST0== tnorm0: %f r0t0: %f tnorm1: %f r0t1: %f tnorm2: %f r0t2: %f   mag: %f\n", tnorm0, r0t0, tnorm1, r0t1, tnorm2, r0t2, sqrt( ((r0[0]-t0[0])*(r0[0]-t0[0])) + ((r0[1]-t0[1])*(r0[1]-t0[1])) + ((r0[2]-t0[2])*(r0[2]-t0[2])) ));
		printf("TEST1== tnorm0: %f r1t0: %f tnorm1: %f r1t1: %f tnorm2: %f r1t2: %f   mag: %f\n", tnorm0, r1t0, tnorm1, r1t1, tnorm2, r1t2, sqrt( ((r1[0]-t0[0])*(r1[0]-t0[0])) + ((r1[1]-t0[1])*(r1[1]-t0[1])) + ((r1[2]-t0[2])*(r1[2]-t0[2])) ));
		exit(1);
	}
	//printf("TRIANGLE CASE: %i  test0: %f test1: %f\n", triangleCase, test0, test1);	///debug instrumentation - DO NOT DELETE
             
	////////return NULL for any cases where the segment is clearly off plane (i.e. no intersection)
	if( triangleCase==0 || triangleCase== 8 ){
		return NULL;
	}

	///Deal with the case where the segment is totally coplanar
	if(triangleCase==4){
		result = intersect_Triangle_Segments(t0, t1, t2, r0, r1);
		return result;
	}
	
	///////////////////////////////////////////////////////////////////////
	///second, identify the point of intersection on the plane.
	double planePt0, planePt1, planePt2;

	///in the cases where there is exactly one point in the plane, this is well defined.
	if( triangleCase==3 || triangleCase==5 ){  ///in this case, r0 is in the plane
		planePt0 = r0[0];		planePt1 = r0[1];		planePt2 = r0[2];
	}
	
	if( triangleCase==1 || triangleCase==7 ){  ///in this case, r1 is in the plane
		planePt0 = r1[0];		planePt1 = r1[1];		planePt2 = r1[2];
	}

	///The most common case is when the points are on either side of the triangle.
	if( triangleCase==2 || triangleCase==6 ){
		///
		///treat the segment as a ray.  We will identify l, the paramaterization variable along this ray,
		///where this ray hits the plane.  We will say the ray hits the plane at point p.
		///
		///This computed by solving two simultaneous equations:
		///Equation 1:      
		///   tnorm dot (p-t0) = 0     "the dot product of the normal and an in-plane vector is zero."
		///
		///Equation 2:
		///   p = r0 + l(r1-r0)        "the point at intersection is along the line of parametrization"
		///
		///Substituting the right side of equation 2 for p in equation 1, gives
		///   tnorm dot ( (r0 + l(r1-r0)) - t0 ) = 0
		///
		///Or, more simply:
		///   tnorm dot( r0 + l(r1-r0) ) = tnorm dot t0
		///
		///Solving for l, this gives
		///
		///          tnorm dot (t0-r0)
		///   l = -----------------------
		///          tnorm dot (r1-r0)
		///

		double topValue = (tnorm0*(t0[0]-r0[0])) + (tnorm1*(t0[1]-r0[1])) + (tnorm2*(t0[2]-r0[2]));
		double botValue = (tnorm0*(r1[0]-r0[0])) + (tnorm1*(r1[1]-r0[1])) + (tnorm2*(r1[2]-r0[2]));
		double l = ( topValue ) / ( botValue );

	//	if( fabs(topValue) < SMALL_NUM ){ topValue = 0; printf("topValue is TINY [%f] (triangle area: %f)\n", topValue, check_area); }
	//	if( fabs(botValue) < SMALL_NUM ){ botValue = 0; printf("botValue is TINY [%f] (triangle area: %f)\n", botValue, check_area); }
	//	if( fabs(l) < SMALL_NUM ){ l = 0; printf("l = [%f] top: [%f] bot: [%f] (triangle area: %f\n", l, topValue, botValue, check_area); }
		
//		///eliminate issues with small numbers
		if(fabs(l) < SMALL_NUM){ test0 = 0.0; }
		
		///if l<0 or l>1, then the parametrization is outside the segment, and they do not touch.
		if( (l>=0) && (l<=1) ){
			///we solve the parametrization, if it is within [0,1]
			///simply plug in l, now that we have it.
			planePt0 = r0[0] + l*(r1[0]-r0[0]);
			planePt1 = r0[1] + l*(r1[1]-r0[1]);
			planePt2 = r0[2] + l*(r1[2]-r0[2]);
		}
		else{
			return NULL;
		}
	}

	///////////////////////////////////////////////////////////////////////
	///third identify if the point is inside the triangle.
	///
	///We will identify if the point is inside the triangle by checking half-planes.
	///If the point is on the correct side of each of the triangle edges, then it is inside.
	///
	///tnorm is (t1-t0) CROSS (t2-t0).  Thus, to determine if the point is
	///inside the triangle, 
	///
	///i0 = tnorm CROSS (t0-t1) gives an in-plane vector pointing inwards of the triangle.
	/// (p-t0) DOT (VECTOR) is positive if the point is on the correct side, negative otherwise.
	///
	///Repeat this process for all edges.
	///
	double hp0, hp1, hp2;
	double tmp0, tmp1, tmp2;

	///crossproduct with first half plane
	tmp0 = ((t0[1]-t1[1])*tnorm2) - ((t0[2]-t1[2])*tnorm1);
	tmp1 = -( ( (t0[0]-t1[0])*tnorm2 ) - ( (t0[2]-t1[2])*tnorm0 ) );
	tmp2 = ((t0[0]-t1[0])*tnorm1) - ((t0[1]-t1[1])*tnorm0);
	hp0 = (tmp0*(planePt0-t0[0])) + (tmp1*(planePt1-t0[1])) + (tmp2*(planePt2-t0[2]));

	///crossproduct with second half plane
	tmp0 = ((t1[1]-t2[1])*tnorm2) - ((t1[2]-t2[2])*tnorm1);
	tmp1 = -( ( (t1[0]-t2[0])*tnorm2 ) - ( (t1[2]-t2[2])*tnorm0 ) );
	tmp2 = ((t1[0]-t2[0])*tnorm1) - ((t1[1]-t2[1])*tnorm0);
	hp1 = (tmp0*(planePt0-t1[0])) + (tmp1*(planePt1-t1[1])) + (tmp2*(planePt2-t1[2]));

	///crossproduct with third half plane
	tmp0 = ((t2[1]-t0[1])*tnorm2) - ((t2[2]-t0[2])*tnorm1);
	tmp1 = -( ( (t2[0]-t0[0])*tnorm2 ) - ( (t2[2]-t0[2])*tnorm0 ) );
	tmp2 = ((t2[0]-t0[0])*tnorm1) - ((t2[1]-t0[1])*tnorm0);
	hp2 = (tmp0*(planePt0-t2[0])) + (tmp1*(planePt1-t2[1])) + (tmp2*(planePt2-t2[2]));

	///fix small-value precision issues
	if(fabs(hp0) < SMALL_NUM){ hp0 = 0.0; }
	if(fabs(hp1) < SMALL_NUM){ hp1 = 0.0; }
	if(fabs(hp2) < SMALL_NUM){ hp2 = 0.0; }

	////this might be a numerical issue with hitting the edge (hp0=0, for example)
	if( (hp0>=0) && (hp1>=0) && (hp2>=0) ){
		result = new double[3];
		result[0] = planePt0;		result[1] = planePt1;		result[2] = planePt2;
	}
	else{
		///if not all cases are true, then we are not inside the triangle.  clear all and return NULL.
		result = NULL;		
	}

	return result;
}








///////////////////////////////////////////////////////////////////////////////////
// segment-Sphere Intersection
///////////////////////////////////////////////////////////////////////////////////
double * intersect_Seg_Sphere( double * r0, double * r1, double * c, double r)
{
	double distance1 = vectorSize(r0[0]-c[0], r0[1]-c[1], r0[2]-c[2]);
	double distance2 = vectorSize(r1[0]-c[0], r1[1]-c[1], r1[2]-c[2]);

	///eliminate anything that doesnt intersect.	
	if( ((distance1>r)&&(distance2>r)) || ((distance1<r)&&(distance2<r)) ){
		return NULL;
	}
	///begin with a segment from r0 to r1
	///
	///The points p on the line are described as p = r0 + u*(r1-r0), or, more specifically,
	///
	///p[0] = r0[0] + u*(r1[0]-r0[0])
	///p[1] = r0[1] + u*(r1[1]-r0[1])
	///p[2] = r0[2] + u*(r1[2]-r0[2])
	///
	///The sphere centered at c with radius r is described as 
	///
	/// (x-c[0])^2 + (y-c[1])^2 + (z - c[2])^2 = r^2
	///
	///If x,y,z is the point of intersection with the line, then x=p[0], y=p[1], z=p[2], and 
	///we can substitute the expressions above into the sphere equation.  If we evaluate the
	///squared quantities ( (x-c[0])^2, etc ) and separate out the terms into terms with
	///u^2, terms with u, and constant terms, we get a quadratic equation of the form
	///
	/// a*(u^2) + b*u + c = 0
	///
	///where 
	///
	/// a = (r1[0]-r0[0])^2 + (r1[1]-r0[1])^2 + (r1[2]-r0[2])^2
	/// b = 2*( ((r1[0]-r0[0])*(r[0]-c[0])) + ((r1[1]-r0[1])*(r[1]-c[1])) + ((r1[2]-r0[2])*(r[2]-c[2])) )
	/// c = c[0]^2 + c[1]^2 + c[2]^2 + r0[0]^2 + r0[1]^2 + r0[2]^2 + 2( (c[0]*r0[0]) + (c[1]*r0[1]) + (c[2]*r0[2]) ) - r^2
	///
	///  Since this is quadratic, its solutions are
	///
	///       -b +-  sqrt( b^2 - 4*a*c )
	///  u = -----------------------------
	///                    2a
	///
	///  Which has a descriminant (b^2 - 4*a*c).  
	///  CASE 1: If the discriminant is less than zero, the segment does not intersect the sphere
	///  CASE 2: If the discriminant is equal to zero, then the segment intersects the sphere at u = (-b)/(2*a)
	///  CASE 3: If the discriminant is greater than zero, then the line intersects the sphere at two places.
	
	///begin by computing the components of the quadratic:
	double quad_a, quad_b, quad_c;
	
	quad_a = ((r1[0]-r0[0])*(r1[0]-r0[0])) + ((r1[1]-r0[1])*(r1[1]-r0[1])) + ((r1[2]-r0[2])*(r1[2]-r0[2]));
	
	quad_b = 2*( ((r1[0]-r0[0])*(r0[0]-c[0])) + ((r1[1]-r0[1])*(r0[1]-c[1])) + ((r1[2]-r0[2])*(r0[2]-c[2])) );
	
	quad_c = (c[0]*c[0]) + (c[1]*c[1]) + (c[2]*c[2]) + (r0[0]*r0[0]) + (r0[1]*r0[1]) + (r0[2]*r0[2]) 
		- (2*( (c[0]*r0[0]) + (c[1]*r0[1]) + (c[2]*r0[2]) )) - (r*r);

	///THis is the quadratic discriminant
	double disc = (quad_b*quad_b) - (4*quad_a*quad_c);
	
	///this is the result
	double * p = NULL;

	//////////////////////////////////////////////////////////////////////////////////////////////
	///case 1 - the line defined by r0 r1 would never touch the sphere
	//////////////////////////////////////////////////////////////////////////////////////////////
	if(disc < -.00001){ 
	//	printf("CASE 1 %f\n", disc);
		p = NULL; 
	}
	//////////////////////////////////////////////////////////////////////////////////////////////
	///case 2 - the line defined by r0 r1 touches the sphere at only one point.  We must find
	///out if it touches 
	//////////////////////////////////////////////////////////////////////////////////////////////
	if( (disc<.00001) && (disc>-.00001) ){
	//	printf("CASE 2 %f\n", disc);
		double u = (-quad_b) / (2*quad_a);
		p = new double[3];
		p[0] = r0[0] + u*(r1[0]-r0[0]);
		p[1] = r0[1] + u*(r1[1]-r0[1]);
		p[2] = r0[2] + u*(r1[2]-r0[2]);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////
	///case 3 - theline defined by r0 r1 touches the sphere at two points.  We must find
	/// out of the line touches within the segment r0 r1
	//////////////////////////////////////////////////////////////////////////////////////////////
	if(disc >.00001){ 
	//	printf("CASE 3 %f\n", disc);
		double u1 = (-quad_b - sqrt( (quad_b*quad_b) - (4*quad_a*quad_c) ) )/(2*quad_a);
		double u2 = (-quad_b + sqrt( (quad_b*quad_b) - (4*quad_a*quad_c) ) )/(2*quad_a);
	
		double u;
		if( u1<0 || u1>1 ){ u = u2; }
		if( u2<0 || u2>1 ){ u = u1; }
	
		p = new double[3];
		p[0] = r0[0] + u*(r1[0]-r0[0]);
		p[1] = r0[1] + u*(r1[1]-r0[1]);
		p[2] = r0[2] + u*(r1[2]-r0[2]);
	}
	
	///return results
	return p;
}






///////////////////////////////////////////////////////////////////////////////////
// tests if a segment AB intersects an axis aligned quad, defined by ranges x1, x2, y1, y2, z1, z2.
// returns 0 if no collision, 1 if collision
///////////////////////////////////////////////////////////////////////////////////
int segAAquadCollision(double * A, double *B,
					double x1, double x2, double y1, double y2, double z1, double z2, 
					double * intersect )
{
	double t;			//t, the parameter
	double xt, yt, zt;	//x, y, z, the position on the segment where it breaks the axis aligned plane.
	if(x1 == x2){
		if( (A[0]-x1)*(B[0]-x1) > 0){ return 0; }
		else{
			t = (A[0]-x1) / (A[0]-B[0]);
			yt = (1-t)*A[1] + t*(B[1]);
			if( yt < y1 || yt > y2 ){ return 0; }
			zt = (1-t)*A[2] + t*(B[2]);
			if( zt < z1 || zt > z2 ){ return 0; }
			xt = x1;	// to fill all values for intersection
		}
	}
	else if(y1 == y2){
		if( (A[1]-y1)*(B[1]-y1) > 0){ return 0; }
		else{
			t = (A[1]-y1) / (A[1]-B[1]);
			xt = (1-t)*A[0] + t*(B[0]);
			if( xt < x1 || xt > x2 ){ return 0; }
			zt = (1-t)*A[2] + t*(B[2]);
			if( zt < z1 || zt > z2 ){ return 0; }
			yt = y1;	// to fill all values for intersection
		}
	}
	else if(z1 == z2){
		if( (A[2]-z1)*(B[2]-z1) > 0){ return 0; }
		else{
			t = (A[2]-z1) / (A[2]-B[2]);
			xt = (1-t)*A[0] + t*(B[0]);
			if( xt < x1 || xt > x2 ){ return 0; }
			yt = (1-t)*A[1] + t*(B[1]);
			if( yt < y1 || yt > y2 ){ return 0; }
			zt = z1;	// to fill all values for intersection
		}
	}

	if( (x1!=x2) && (y1!=y2) && (z1!=z2) ){
		printf("ERROR [segAAquadCollision()]: Axis aligned quad is not axis aligned! exitting.\n");
		exit(1);
	}

	if( intersect != NULL){
		intersect[0] = xt;
		intersect[1] = yt;
		intersect[2] = zt;
	}

	return 1;
}



///////////////////////////////////////////////////////////////////////////////////
///determine if a triangle intersects the volume of the cube. (tcoords is 9 doubles = 3*3)
// this code uses the following cube layout, with 3D axes to the right.
//      7---6            ^
//      |\  |\       \   |
//      | 4---5       +y |
//      | | | |        \ z+
//      3-|-2 |         \|  
//       \|  \|          0--x+-->    
//        0---1
///////////////////////////////////////////////////////////////////////////////////
bool triCubeIntersect(double * tcoords, double xn, double xp, double yn, double yp, double zn, double zp)
{
	int i = 0;
	
	///this gets tested in two cases.
	bool result = false;
	
	//First, test if any one of the triangle endpoints is in the cube.
	///if any point is within the bounds, the triangle is within the cube.
	for(i = 0; i<3; i++){
		if(  (tcoords[3*i+0] > xn) && (tcoords[3*i+0] < xp) &&
			(tcoords[3*i+1] > yn) && (tcoords[3*i+1] < yp) &&
			(tcoords[3*i+2] > zn) && (tcoords[3*i+2] < zp) )
		{
			return true;
		}
	}
	
	double ** cc = new double*[8];
	cc[0] = new double[3];  cc[0][0] = xn;  cc[0][1] = yn;  cc[0][2] = zn;  //point 0   this code uses the following cube layout
	cc[1] = new double[3];  cc[1][0] = xp;  cc[1][1] = yn;  cc[1][2] = zn;  //point 1         7---6            ^                     
	cc[2] = new double[3];  cc[2][0] = xp;  cc[2][1] = yp;  cc[2][2] = zn;  //point 2         |\  |\       \   |                     
	cc[3] = new double[3];  cc[3][0] = xn;  cc[3][1] = yp;  cc[3][2] = zn;  //point 3         | 4---5       +y |                     
	cc[4] = new double[3];  cc[4][0] = xn;  cc[4][1] = yn;  cc[4][2] = zp;  //point 4         | | | |        \ z+                    
	cc[5] = new double[3];  cc[5][0] = xp;  cc[5][1] = yn;  cc[5][2] = zp;  //point 5         3-|-2 |         \|                     
	cc[6] = new double[3];  cc[6][0] = xp;  cc[6][1] = yp;  cc[6][2] = zp;  //point 6          \|  \|          0--x+-->              
	cc[7] = new double[3];  cc[7][0] = xn;  cc[7][1] = yp;  cc[7][2] = zp;  //point 7           0---1                                
	
	if(!result){ result = triQuadOverlapTest( tcoords, cc[0], cc[3], cc[2], cc[1] ); }
	if(!result){ result = triQuadOverlapTest( tcoords, cc[0], cc[1], cc[5], cc[4] ); }
	if(!result){ result = triQuadOverlapTest( tcoords, cc[1], cc[2], cc[6], cc[5] ); }
	if(!result){ result = triQuadOverlapTest( tcoords, cc[2], cc[3], cc[7], cc[6] ); }
	if(!result){ result = triQuadOverlapTest( tcoords, cc[3], cc[0], cc[4], cc[7] ); }
	if(!result){ result = triQuadOverlapTest( tcoords, cc[4], cc[5], cc[6], cc[7] ); }
	
	//clean up
	for(i = 0; i<8; i++){ delete[]( cc[i]); }
	delete[](cc);
	
	return result;
}











///////////////////////////////////////////////////////////////////////////////////
///wrapper added by Brian Chen
///	note that order of corner enumeration is important!!!!
///	It is assumed that the four corners make a circuit around the quad 
///	i.e. q1, q2, q3, q4 draws the quad, not a "bow-tie"
///////////////////////////////////////////////////////////////////////////////////
bool triQuadOverlapTest( double tcoords[9], double q0[3], double q1[3], double q2[3], double q3[3] )
{
	double * t0 = new double[3]; t0[0] = tcoords[0]; t0[1] = tcoords[1]; t0[2] = tcoords[2];
	double * t1 = new double[3]; t1[0] = tcoords[3]; t1[1] = tcoords[4]; t1[2] = tcoords[5];
	double * t2 = new double[3]; t2[0] = tcoords[6]; t2[1] = tcoords[7]; t2[2] = tcoords[8];

	///test overlaps.  make sure triangles are right.
	int info1 = tri_tri_overlap_test_3d(t0, t1, t2, q0, q1, q2);
	int info2 = tri_tri_overlap_test_3d(t0, t1, t2, q2, q3, q0);

	delete[](t0);delete[](t1);delete[](t2);

	///test if either overlaps
	if( (info1==1) || (info2==1) ){
		return true;
	}
	return false;
}






/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/




