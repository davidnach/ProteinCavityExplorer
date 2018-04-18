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
 *       surface generation for abstract spheres.
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


#include "primativesProcessing.h"

class Primative;
class PrimSphere;
class PrimCube;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////void constructor for Primative
Primative::Primative(){ }

////void destructor for Primative
Primative::~Primative(){ }

////virtual function does nothing
bool Primative::isInside(double * pt){ }

////virtual function does nothing
double * Primative::intersectionPt(double * pt1, double * pt2){ }

////virtual function does nothing
double * Primative::getBounds(){ }

////virtual function does nothing
Primative * Primative::copy(){ }

////virtual function does nothing
SurfaceObject * Primative::genSurface(){ }

////virtual function does nothing
double * Primative::normal(double * pt){ }

////virtual function does nothing
double Primative::distance(double * pt){ }



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////General intersection test.
//
//
bool Primative::intersects( Primative * p )
{
	if (primType == PRIM_TYPE_SPHERE){
		if(p->primType == PRIM_TYPE_SPHERE){
			return sphereIntersectsSphere( (PrimSphere*) this, (PrimSphere*) p);
		}
		if(p->primType == PRIM_TYPE_CUBE){
			return sphereIntersectsCube( (PrimSphere*) this, (PrimCube*) p);
		}
	}
	if (primType == PRIM_TYPE_CUBE){
		if(p->primType == PRIM_TYPE_SPHERE){
			return sphereIntersectsCube( (PrimSphere*) p, (PrimCube*) this);
		}
		if(p->primType == PRIM_TYPE_CUBE){
			return cubeIntersectsCube( (PrimCube*) this, (PrimCube*) p);
		}
	}
	
	printf("ERROR: Primative Intersection Error.  This should never happen!\n");
	return false;	
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////Sphere-sphere specific intersection test
//
//
bool Primative::sphereIntersectsSphere( PrimSphere * s1, PrimSphere * s2)
{
	bool result = false;
	double dist = vectorSize( s1->origin[0] - s2->origin[0], s1->origin[1] - s2->origin[1], s1->origin[2] - s2->origin[2] );
	if( dist <= ( s1->radius + s2->radius ) ){
		result = true;
	}
	
	return result;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////Sphere-cube specific intersection test
//
//
bool Primative::sphereIntersectsCube  ( PrimSphere * s1, PrimCube * c2)
{
	///if outside the bounding cube.
	bool test1 = s1->origin[0] < ((c2->center[0] - c2->size) - s1->radius);
	bool test2 = s1->origin[0] > ((c2->center[0] + c2->size) + s1->radius);
	bool test3 = s1->origin[1] < ((c2->center[1] - c2->size) - s1->radius);
	bool test4 = s1->origin[1] > ((c2->center[1] + c2->size) + s1->radius);
	bool test5 = s1->origin[2] < ((c2->center[2] - c2->size) - s1->radius);
	bool test6 = s1->origin[2] > ((c2->center[2] + c2->size) + s1->radius);
	
	if( test1 || test2 || test3 || test4 || test5 || test6 ){ return false; }
	
	///if inside the input cube.
	int test1a; if( s1->origin[0] <= (c2->center[0] + c2->size) ){ test1a = 1; } else{ test1a = 0; }
	int test2a; if( s1->origin[0] >= (c2->center[0] - c2->size) ){ test2a = 1; } else{ test2a = 0; }
	int test3a; if( s1->origin[1] <= (c2->center[1] + c2->size) ){ test3a = 1; } else{ test3a = 0; }
	int test4a; if( s1->origin[1] >= (c2->center[1] - c2->size) ){ test4a = 1; } else{ test4a = 0; }
	int test5a; if( s1->origin[2] <= (c2->center[2] + c2->size) ){ test5a = 1; } else{ test5a = 0; }
	int test6a; if( s1->origin[2] >= (c2->center[2] - c2->size) ){ test6a = 1; } else{ test6a = 0; }
	
	int sum = test1a + test2a + test3a + test4a + test5a + test6a;
	if( sum == 6 ){ return true; }	///the circle origin is within the cube or on its surface.
	if( sum == 5 ){ return true; }	///the circle origin is between a face of the cube and the bounding cube.
	if( sum == 4 ){				///test if the circle origin is in the cylinder along an edge of the cube.
		double cd1plus, cd1minus, cd2plus, cd2minus;
		if( test1a + test2a == 2 ){	//not in the X direction.
			cd1plus  = c2->center[1]+c2->size;	cd1minus = c2->center[1]-c2->size;
			cd2plus  = c2->center[2]+c2->size;	cd2minus = c2->center[2]-c2->size;
		}
		if( test3a + test4a == 2 ){
			cd1plus  = c2->center[0]+c2->size;	cd1minus = c2->center[0]-c2->size;
			cd2plus  = c2->center[2]+c2->size;	cd2minus = c2->center[2]-c2->size;
		}
		if( test5a + test6a == 2 ){
			cd1plus  = c2->center[0]+c2->size;	cd1minus = c2->center[0]-c2->size;
			cd2plus  = c2->center[1]+c2->size;	cd2minus = c2->center[1]-c2->size;
		}
		double d1 = vectorSize2d( (cd1plus)-s1->origin[1]  , (cd2plus)-s1->origin[2]  ); if(d1<=s1->radius){return true;}
		double d2 = vectorSize2d( (cd1plus)-s1->origin[1]  , (cd2minus)-s1->origin[2] ); if(d2<=s1->radius){return true;}
		double d3 = vectorSize2d( (cd1minus)-s1->origin[1] , (cd2plus)-s1->origin[2]  ); if(d3<=s1->radius){return true;}
		double d4 = vectorSize2d( (cd1minus)-s1->origin[1] , (cd2minus)-s1->origin[2] ); if(d4<=s1->radius){return true;}
	}

	if( sum == 3 ){				///test if the circle origin is in the sphere corner at the edge of the cube.
		double as0 = s1->origin[0];	double as1 = s1->origin[1];	double as2 = s1->origin[2];
		double cd0p = c2->center[0]+c2->size; double cd0m = c2->center[0]-c2->size; 
		double cd1p = c2->center[1]+c2->size; double cd1m = c2->center[1]-c2->size; 
		double cd2p = c2->center[2]+c2->size; double cd2m = c2->center[2]-c2->size; 
		double d;
		d = vectorSize(cd0p-as0, cd1p-as1, cd2p-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0p-as0, cd1p-as1, cd2m-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0p-as0, cd1m-as1, cd2p-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0p-as0, cd1m-as1, cd2m-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0m-as0, cd1p-as1, cd2p-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0m-as0, cd1p-as1, cd2m-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0m-as0, cd1m-as1, cd2p-as2); if(d<=s1->radius){ return true; }
		d = vectorSize(cd0m-as0, cd1m-as1, cd2m-as2); if(d<=s1->radius){ return true; }
		
	}

	return false;

}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////cube-cube intersectionTest
//
//
bool Primative::cubeIntersectsCube    ( PrimCube   * c1, PrimCube * c2)
{
	if( ( c1->center[0]+c1->size < c2->center[0]-c2->size) || ( c1->center[0]-c1->size > c2->center[0]+c2->size) ) { return false;}
	if( ( c1->center[1]+c1->size < c2->center[1]-c2->size) || ( c1->center[1]-c1->size > c2->center[1]+c2->size) ) { return false;}
	if( ( c1->center[2]+c1->size < c2->center[2]-c2->size) || ( c1->center[2]-c1->size > c2->center[2]+c2->size) ) { return false;}
	
	return true;	
}	



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//Primitive Sphere constructor
//
//
PrimSphere::PrimSphere(double * o, double r)
{
	origin = new double[3];
	origin[0] = o[0];
	origin[1] = o[1];
	origin[2] = o[2];

	radius = r;
	primType = PRIM_TYPE_SPHERE;
}	

/////////////////////////////////////////////////////////////////////////////
////destructor for primative sphere	
PrimSphere::~PrimSphere()
{
	delete[](origin);
}

/////////////////////////////////////////////////////////////////////////////
////Returns the dimensions of the bounding box
double * PrimSphere::getBounds()
{
	double * result = new double[6];
	result[0] = origin[0]-radius;
	result[1] = origin[0]+radius;
	result[2] = origin[1]-radius;
	result[3] = origin[1]+radius;
	result[4] = origin[2]-radius;
	result[5] = origin[2]+radius;

	return result;
}


	
/////////////////////////////////////////////////////////////////////////////
////True if pt is inside the sphere, false otherwise.
bool PrimSphere::isInside(double * pt)
{
	bool result = false;
	double dist = vectorSize( pt[0]-origin[0], pt[1]-origin[1], pt[2]-origin[2] );
	if( dist <= radius ){
		result = true;
	}
	
	return result;
}


/////////////////////////////////////////////////////////////////////////////
////Point of intersection of the segment with the sphere.
///it is understood that pt1 is always inside
double * PrimSphere::intersectionPt(double * pt1, double * pt2)
{
	double * result = intersect_Seg_Sphere( pt1, pt2, origin, radius);
	return result;
}


/////////////////////////////////////////////////////////////////////////////
/////Create a copy of this primSphere
PrimSphere * PrimSphere::copy()
{
	PrimSphere * newSphere = new PrimSphere(origin, radius);
	return newSphere;
}

/////////////////////////////////////////////////////////////////////////////
////creates a procedural surface for this sphere.
///currently this does nothing.  Dont call it.
SurfaceObject * PrimSphere::genSurface()
{
	printf("ERROR: PrimSphere currently does not generate procedural surfaces.\n");
	printf("ERROR: returning NULL and exitting.\n");
	exit(1);
	return NULL;
	
}


/////////////////////////////////////////////////////////////////////////////
////returns the normal 
///currently this does nothing.  Dont call it.
double * PrimSphere::normal(double * pt)
{
	///nothing here
	double * result = new double[3];
	result[0] = pt[0]-origin[0];
	result[1] = pt[1]-origin[1];
	result[2] = pt[2]-origin[2];
	
	return result;
}


/////////////////////////////////////////////////////////////////////////////
////returns the normal 
///
double PrimSphere::distance(double * pt)
{
	double dist = vectorSize( pt[0]-origin[0], pt[1]-origin[1], pt[2]-origin[2] );

	///assuming dist < radius
	double result = radius-dist;
	///else
	if(dist > radius){ result = dist-radius; }

	return result;
	
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//  Primative Cube Constructor
//	double * center;	///double[3], x, y, z - center of cube
//	double size;		///one half the sidelength of the cube.
PrimCube::PrimCube(double * c, double s)
{
	center = new double[3];
	center[0] = c[0];
	center[1] = c[1];
	center[2] = c[2];
	size = s;
	primType = PRIM_TYPE_CUBE;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//  Primative Cube destructor
PrimCube::~PrimCube()
{
	delete[](center);
}


/////////////////////////////////////////////////////////////////////////////
////Returns the dimensions of the bounding box
double * PrimCube::getBounds()
{
	double * result = new double[6];
	result[0] = center[0]-size;
	result[1] = center[0]+size;
	result[2] = center[1]-size;
	result[3] = center[1]+size;
	result[4] = center[2]-size;
	result[5] = center[2]+size;

	return result;
}
	
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//  Primative Cube inside test
bool PrimCube::isInside(double * pt)
{
//	printf("cube isInside()\n");
	
	if( (pt[0] < center[0]-size) || (pt[0] > center[0]+size) ){ return false; }
	if( (pt[1] < center[1]-size) || (pt[1] > center[1]+size) ){ return false; }
	if( (pt[2] < center[2]-size) || (pt[2] > center[2]+size) ){ return false; }
	
	return true;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//  Primative Cube intersection test
///
/// NOTE: THIS IS NOT A GENERAL SEGMENT-CUBE INTERSECTION TEST.  IT RETURNS EXACTLY ONE POINT
/// BECAUSE IT IS ASSUMING THAT ONE POINT IS INSIDE AND ONE OUTSIDE
///
/// IT IS ALSO ASSUMING THAT THIS IS AN AXIS ALIGNED CUBE, WHICH IS THIS TYPE.
/// IT IS ALSO ASSUMING THAT ALL SEGMENTS POSED TO IT ARE AXIS ALIGNED CUBES, SINCE WE ONLY USE IT FOR MC
double * PrimCube::intersectionPt(double * pt1, double * pt2)
{
	if( !isInside(pt1) && !isInside(pt2) ){
		return NULL;
	}
	
	double * result = new double[3];
	int dimState = -1;
	bool AAcount = 0;
	if(pt1[0] == pt2[0]){ result[0] = pt1[0]; }	else{ dimState = 0; AAcount++; }
	if(pt1[1] == pt2[1]){ result[1] = pt1[1]; }	else{ dimState = 1; AAcount++; }
	if(pt1[2] == pt2[2]){ result[2] = pt1[2]; }	else{ dimState = 2; AAcount++; }
	if(AAcount > 1){
		printf("ERROR: PrimCube::intersectionPt() was given a non axis aligned segment! exitting!\n");
		exit(1);
	}

	double loVal = pt1[dimState];
	if(pt1[dimState] > pt2[dimState]){ loVal = pt2[dimState]; }
	double hiVal = pt1[dimState];
	if(pt1[dimState] < pt2[dimState]){ hiVal = pt2[dimState]; }

	if( (loVal <= center[dimState]-size) && (hiVal >= center[dimState]-size) ){ result[dimState] = center[dimState]-size; }
	if( (loVal <= center[dimState]+size) && (hiVal >= center[dimState]+size) ){ result[dimState] = center[dimState]+size; }

	return result;
}


/////////////////////////////////////////////////////////////////////////////
/////Create a copy of this primSphere
PrimCube * PrimCube::copy()
{
	PrimCube * newCube = new PrimCube(center, size);
	return newCube;
}

/////////////////////////////////////////////////////////////////////////////
////creates a procedural surface for this sphere.
///currently this does nothing.  Dont call it.
SurfaceObject * PrimCube::genSurface()
{
	int surf_numPts = 8;
	double * surf_pts = new double[8*3];
	double * surf_norms = new double[8*3];
	int surf_numTris = 12;
	int * surf_tris = new int[12*3];
	double * surf_triangleNorms = new double[12*3];

	surf_pts[3*0+0]=center[0]-size; surf_pts[3*0+1]=center[1]-size; surf_pts[3*0+2]=center[2]-size; 
	surf_pts[3*1+0]=center[0]-size; surf_pts[3*1+1]=center[1]-size; surf_pts[3*1+2]=center[2]+size; 
	surf_pts[3*2+0]=center[0]-size; surf_pts[3*2+1]=center[1]+size; surf_pts[3*2+2]=center[2]-size; 
	surf_pts[3*3+0]=center[0]-size; surf_pts[3*3+1]=center[1]+size; surf_pts[3*3+2]=center[2]+size; 
	surf_pts[3*4+0]=center[0]+size; surf_pts[3*4+1]=center[1]-size; surf_pts[3*4+2]=center[2]-size; 
	surf_pts[3*5+0]=center[0]+size; surf_pts[3*5+1]=center[1]-size; surf_pts[3*5+2]=center[2]+size; 
	surf_pts[3*6+0]=center[0]+size; surf_pts[3*6+1]=center[1]+size; surf_pts[3*6+2]=center[2]-size; 
	surf_pts[3*7+0]=center[0]+size; surf_pts[3*7+1]=center[1]+size; surf_pts[3*7+2]=center[2]+size; 

	surf_norms[3*0+0]= -1.0 ; surf_norms[3*0+1]= -1.0 ; surf_norms[3*0+2]= -1.0 ; ///diagonal norms
	surf_norms[3*1+0]= -1.0 ; surf_norms[3*1+1]= -1.0 ; surf_norms[3*1+2]= +1.0 ; 
	surf_norms[3*2+0]= -1.0 ; surf_norms[3*2+1]= +1.0 ; surf_norms[3*2+2]= -1.0 ; 
	surf_norms[3*3+0]= -1.0 ; surf_norms[3*3+1]= +1.0 ; surf_norms[3*3+2]= +1.0 ; 
	surf_norms[3*4+0]= +1.0 ; surf_norms[3*4+1]= -1.0 ; surf_norms[3*4+2]= -1.0 ; 
	surf_norms[3*5+0]= +1.0 ; surf_norms[3*5+1]= -1.0 ; surf_norms[3*5+2]= +1.0 ; 
	surf_norms[3*6+0]= +1.0 ; surf_norms[3*6+1]= +1.0 ; surf_norms[3*6+2]= -1.0 ; 
	surf_norms[3*7+0]= +1.0 ; surf_norms[3*7+1]= +1.0 ; surf_norms[3*7+2]= +1.0 ; 

	surf_tris[3*0+0] =  0; surf_tris[3*0+1] =  2; surf_tris[3*0+2] =   6;   //floor   0, 2, 6  
	surf_tris[3*1+0] =  0; surf_tris[3*1+1] =  6; surf_tris[3*1+2] =   4;   //floor   0, 6, 4  
	surf_tris[3*2+0] =  1; surf_tris[3*2+1] =  5; surf_tris[3*2+2] =   7;   //cieling 1, 5, 7  
	surf_tris[3*3+0] =  1; surf_tris[3*3+1] =  7; surf_tris[3*3+2] =   3;   //cieling 1, 7, 3  
	surf_tris[3*4+0] =  0; surf_tris[3*4+1] =  1; surf_tris[3*4+2] =   3;   //left    0, 1, 3, 
	surf_tris[3*5+0] =  0; surf_tris[3*5+1] =  3; surf_tris[3*5+2] =   2;   //left    0, 3, 2, 
	surf_tris[3*6+0] =  4; surf_tris[3*6+1] =  6; surf_tris[3*6+2] =   7;   //right   4, 6, 7, 
	surf_tris[3*7+0] =  4; surf_tris[3*7+1] =  7; surf_tris[3*7+2] =   5;   //right   4, 7, 5, 
	surf_tris[3*8+0] =  0; surf_tris[3*8+1] =  5; surf_tris[3*8+2] =   1;   //front   0, 5, 1, 
	surf_tris[3*9+0] =  0; surf_tris[3*9+1] =  4; surf_tris[3*9+2] =   5;   //front   0, 4, 5, 
	surf_tris[3*10+0] = 2; surf_tris[3*10+1] = 3; surf_tris[3*10+2] =  7;   //back    2, 3, 7, 
	surf_tris[3*11+0] = 2; surf_tris[3*11+1] = 7; surf_tris[3*11+2] =  6;   //back    2, 7, 6, 

	surf_triangleNorms[3*0+0] =  0.0; surf_triangleNorms[3*0+1] =  0.0; surf_triangleNorms[3*0+2] =  -1.0;    //floor   0, 2, 6  
	surf_triangleNorms[3*1+0] =  0.0; surf_triangleNorms[3*1+1] =  0.0; surf_triangleNorms[3*1+2] =  -1.0;    //floor   0, 6, 4  
	surf_triangleNorms[3*2+0] =  0.0; surf_triangleNorms[3*2+1] =  0.0; surf_triangleNorms[3*2+2] =   1.0;    //cieling 1, 5, 7  
	surf_triangleNorms[3*3+0] =  0.0; surf_triangleNorms[3*3+1] =  0.0; surf_triangleNorms[3*3+2] =   1.0;    //cieling 1, 7, 3  
	surf_triangleNorms[3*4+0] = -1.0; surf_triangleNorms[3*4+1] =  0.0; surf_triangleNorms[3*4+2] =   0.0;    //left    0, 1, 3, 
	surf_triangleNorms[3*5+0] = -1.0; surf_triangleNorms[3*5+1] =  0.0; surf_triangleNorms[3*5+2] =   0.0;    //left    0, 3, 2, 
	surf_triangleNorms[3*6+0] =  1.0; surf_triangleNorms[3*6+1] =  0.0; surf_triangleNorms[3*6+2] =   0.0;    //right   4, 6, 7, 
	surf_triangleNorms[3*7+0] =  1.0; surf_triangleNorms[3*7+1] =  0.0; surf_triangleNorms[3*7+2] =   0.0;    //right   4, 7, 5, 
	surf_triangleNorms[3*8+0] =  0.0; surf_triangleNorms[3*8+1] = -1.0; surf_triangleNorms[3*8+2] =   0.0;    //front   0, 5, 1, 
	surf_triangleNorms[3*9+0] =  0.0; surf_triangleNorms[3*9+1] = -1.0; surf_triangleNorms[3*9+2] =   0.0;    //front   0, 4, 5, 
	surf_triangleNorms[3*10+0] = 0.0; surf_triangleNorms[3*10+1] = 1.0; surf_triangleNorms[3*10+2] =  0.0;    //back    2, 3, 7, 
	surf_triangleNorms[3*11+0] = 0.0; surf_triangleNorms[3*11+1] = 1.0; surf_triangleNorms[3*11+2] =  0.0;    //back    2, 7, 6, 

//	SurfaceObject::SurfaceObject(int numPts, double * pts, double * norms, int numTris, int * tris, double * triangleNorms)
	SurfaceObject * result = new SurfaceObject(surf_numPts, surf_pts, surf_norms, surf_numTris, surf_tris, surf_triangleNorms);

	delete[](surf_pts);
	delete[](surf_norms);
	delete[](surf_tris);
	delete[](surf_triangleNorms);

	return result;

}


/////////////////////////////////////////////////////////////////////////////
////returns the normal 
double * PrimCube::normal(double * pt)
{
	double * p0 = new double[3]; p0[0]=center[0]-size; p0[1]=center[1]-size; p0[2]=center[2]-size; 
	double * p1 = new double[3]; p1[0]=center[0]-size; p1[1]=center[1]-size; p1[2]=center[2]+size; 
	double * p2 = new double[3]; p2[0]=center[0]-size; p2[1]=center[1]+size; p2[2]=center[2]-size; 
	double * p3 = new double[3]; p3[0]=center[0]-size; p3[1]=center[1]+size; p3[2]=center[2]+size; 
	double * p4 = new double[3]; p4[0]=center[0]+size; p4[1]=center[1]-size; p4[2]=center[2]-size; 
	double * p5 = new double[3]; p5[0]=center[0]+size; p5[1]=center[1]-size; p5[2]=center[2]+size; 
	double * p6 = new double[3]; p6[0]=center[0]+size; p6[1]=center[1]+size; p6[2]=center[2]-size; 
	double * p7 = new double[3]; p7[0]=center[0]+size; p7[1]=center[1]+size; p7[2]=center[2]+size; 

	double minDistance = HUGE_VAL;
	double tempDist;
	double * result = new double[3];
	
	tempDist = distanceTo(pt, p0, p2, p6); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] =  0.0; result[2] = -1.0; }
	tempDist = distanceTo(pt, p0, p6, p4); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] =  0.0; result[2] = -1.0; }
	tempDist = distanceTo(pt, p1, p5, p7); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] =  0.0; result[2] =  1.0; }
	tempDist = distanceTo(pt, p1, p7, p3); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] =  0.0; result[2] =  1.0; }
	tempDist = distanceTo(pt, p0, p1, p3); if(tempDist < minDistance){ minDistance = tempDist; result[0] = -1.0; result[1] =  0.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p0, p3, p2); if(tempDist < minDistance){ minDistance = tempDist; result[0] = -1.0; result[1] =  0.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p4, p6, p7); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  1.0; result[1] =  0.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p4, p7, p5); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  1.0; result[1] =  0.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p0, p5, p1); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] = -1.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p0, p4, p5); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] = -1.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p2, p3, p7); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] =  1.0; result[2] =  0.0; }
	tempDist = distanceTo(pt, p2, p7, p6); if(tempDist < minDistance){ minDistance = tempDist; result[0] =  0.0; result[1] =  1.0; result[2] =  0.0; }

	delete[](p0);	delete[](p1);	delete[](p2);	delete[](p3);	delete[](p4);	delete[](p5);	delete[](p6);	delete[](p7);

	return result;
}



/////////////////////////////////////////////////////////////////////////////
////returns the normal 
///
double PrimCube::distance(double * pt)
{
	double * p0 = new double[3]; p0[0]=center[0]-size; p0[1]=center[1]-size; p0[2]=center[2]-size; 
	double * p1 = new double[3]; p1[0]=center[0]-size; p1[1]=center[1]-size; p1[2]=center[2]+size; 
	double * p2 = new double[3]; p2[0]=center[0]-size; p2[1]=center[1]+size; p2[2]=center[2]-size; 
	double * p3 = new double[3]; p3[0]=center[0]-size; p3[1]=center[1]+size; p3[2]=center[2]+size; 
	double * p4 = new double[3]; p4[0]=center[0]+size; p4[1]=center[1]-size; p4[2]=center[2]-size; 
	double * p5 = new double[3]; p5[0]=center[0]+size; p5[1]=center[1]-size; p5[2]=center[2]+size; 
	double * p6 = new double[3]; p6[0]=center[0]+size; p6[1]=center[1]+size; p6[2]=center[2]-size; 
	double * p7 = new double[3]; p7[0]=center[0]+size; p7[1]=center[1]+size; p7[2]=center[2]+size; 

	double minDistance = HUGE_VAL;
	double tempDist;
	
	tempDist = distanceTo(pt, p0, p2, p6); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p0, p6, p4); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p1, p5, p7); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p1, p7, p3); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p0, p1, p3); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p0, p3, p2); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p4, p6, p7); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p4, p7, p5); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p0, p5, p1); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p0, p4, p5); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p2, p3, p7); if(tempDist < minDistance){ minDistance = tempDist; }
	tempDist = distanceTo(pt, p2, p7, p6); if(tempDist < minDistance){ minDistance = tempDist; }

	delete[](p0);	delete[](p1);	delete[](p2);	delete[](p3);	delete[](p4);	delete[](p5);	delete[](p6);	delete[](p7);

	return minDistance;
}



	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///Constructor for primativesSet
//	set_t objs;
primativesSet::primativesSet()
{
	objs = alloc_set(SP_MAP);
	isInsideOut = false;
	hash = NULL;
}
	

primativesSet::primativesSet(char * filename)
{
	isInsideOut = false;
	objs = alloc_set(SP_MAP);

	FILE * file = fopen(filename, "r");
	if(file == NULL){ printf("Volumetric Primatives File [%s] Not Found.\n", filename); exit(1); }

	///Set up parsing buffers.
	char * line = new char[1000];
	char * ptr;
	
	///get the points out of the file
	while( fgets(line, NUMCOLUMNS, file) != NULL){
		if(line[0] == '#'){continue;}
		if(strlen(line) == 0){continue;}
		else{
			ptr = strtok(line, " \r\t\n");
			if(ptr == NULL){ continue; }
			if(strcmp(ptr, "SPHERE")==0){
				double * c = new double[3];
				double * r = new double[1];
				c[0] = atof(strtok(NULL, " \r\t\n"));
				c[1] = atof(strtok(NULL, " \r\t\n"));
				c[2] = atof(strtok(NULL, " \r\t\n"));
				r[0] = atof(strtok(NULL, " \r\t\n"));
				printf("Found Sphere: %f %f %f  Radius: %f\n", c[0], c[1], c[2], r[0]);
				
				PrimSphere * tempSphere = new PrimSphere(c, r[0]);
				addObj( (Primative *) tempSphere);
				delete[](c);
				delete[](r);
				continue;
			}
			
			if(strcmp(ptr, "CUBE")==0){
				double * c = new double[3];
				double * r = new double[1];
				c[0] = atof(strtok(NULL, " \r\t\n"));
				c[1] = atof(strtok(NULL, " \r\t\n"));
				c[2] = atof(strtok(NULL, " \r\t\n"));
				r[0] = atof(strtok(NULL, " \r\t\n"));
				printf("Found Cube: %f %f %f  halfSide: %f\n", c[0], c[1], c[2], r[0]);
				
				PrimCube * tempCube = new PrimCube(c, r[0]);
				addObj( (Primative *) tempCube);
				delete[](c);
				delete[](r);
				continue;
			}
		}
	}	

	delete[](line);
	
	hash = NULL;

}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///Destructor for primativesSet
//	set_t objs;
primativesSet::~primativesSet()
{
	int i = 0;
	for(i = 0; i<size_set(objs); i++){
		Primative * temp = (Primative *) mapsto_set(objs, i);
		delete(temp);
	}
	free_set(objs);
}
	
int primativesSet::size()
{
	return size_set(objs);
}

///returns a pointer to the obj desired; NULL if not in range.
Primative * primativesSet::getPrim(int ind)
{
	if( (ind < 0) || (ind >= size()) ){ return NULL; }

	Primative * temp = (Primative *) mapsto_set(objs, ind);
	return temp;
}


void primativesSet::addObj(Primative * obj)
{
	objs = associate_set(objs, size_set(objs), obj);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
////This version of isInside doesnt go through all primatives in the set, under the assumption that there
////are many many primatives.  Instead, if uses a LatticeHash as a member of this object to hash the locations
////and permits a more specific search.
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//dispatch function - calls the appropriate Function
bool primativesSet::isInside(double * pt)
{
	if( hash == NULL ){
		return isInsideLinear(pt);
	}
	else{
		return isInsideHashed(pt);
	}
}
	
//default; checks all primatives for every point.
bool primativesSet::isInsideLinear(double * pt)
{
	bool result = false;
	int i = 0;
	for(i = 0; i<size_set(objs); i++){
		Primative * temp = (Primative *) mapsto_set(objs, i);
		if( temp->isInside(pt) ){
			result = true;
			break;
		}
	}
	
	if(isInsideOut){
		return !result;
	}
	return result;	
}


//uses a hash table to check only nearby primatives.
bool primativesSet::isInsideHashed(double * pt)
{
	bool result = false;

//	printf("OMG WTF BBQ");

	// get objects in the same cube as a point, appends them to list
	set_t listOfPrims = alloc_set(0);
	listOfPrims = hash->getNearbyObjs( pt, listOfPrims );

	int i = 0;
	for(i = 0; i<size_set(listOfPrims); i++){
		Primative * temp = (Primative *) mapsto_set(objs, listOfPrims[i]);
		if( temp->isInside(pt) ){
			result = true;
			break;
		}
	}
	
	free_set(listOfPrims);
	
	if(isInsideOut){
		return !result;
	}
	return result;	
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// This version of intersectionPt uses a LatticeHash to check only the cubes that intersec the
// segment, so that a minimum number of cubes is tested, pointing to a minimum number of possible
// primatives to test.
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//this function dispatches to the other two.
double * primativesSet::intersectionPt(double * pt1, double * pt2)
{
	if( hash == NULL ){
		return intersectionPtLinear(pt1, pt2);
	}
	else{
		return intersectionPtHashed(pt1, pt2);
	}
}


//finds the average intersedction point along a segment with the primatives the segment intersects
double * primativesSet::intersectionPtLinear(double * pt1, double * pt2)
{
	int i = 0;
	double * result = new double[3];
	result[0] = 0;
	result[1] = 0;
	result[2] = 0;
	int counter = 0;
	for(i = 0; i<size_set(objs); i++){
		Primative * temp = (Primative *) mapsto_set(objs, i);
		double * tempResult = temp->intersectionPt(pt1, pt2);
		if(tempResult != NULL){
			result[0] += tempResult[0];
			result[1] += tempResult[1];
			result[2] += tempResult[2];
			delete[](tempResult);
			counter++;
		}
	}
	result[0] /= counter;
	result[1] /= counter;
	result[2] /= counter;
	
	return result;
	
}

//finds the intersedction point of a segment that is inside all the primsets and outside the others.
double * primativesSet::intersectionPtHashed(double * pt1, double * pt2)
{
	// get objects in the same cubes as a segment, appends them to list
	set_t list = alloc_set(0);
	list = hash->getNearbyObjs( pt1, pt2, list );

	int i = 0;
	double * result = new double[3];
	result[0] = 0;
	result[1] = 0;
	result[2] = 0;
	int counter = 0;
	
	for(i = 0; i<size_set(list); i++){
		Primative * temp = (Primative *) mapsto_set(objs, list[i]);
		double * tempResult = temp->intersectionPt(pt1, pt2);
		if(tempResult != NULL){
			result[0] += tempResult[0];
			result[1] += tempResult[1];
			result[2] += tempResult[2];
			delete[](tempResult);
			counter++;
		}
	}
	result[0] /= counter;
	result[1] /= counter;
	result[2] /= counter;
	
	free_set(list);

	return result;
		
	
	
}







double * primativesSet::getPrimSetBounds()
{
	if( size() == 0){
		printf("GETPRIMSETBOUNDS: Empty primatives set! Nothing to do! exitting.\n");
		exit(1);
	}
	
	int i = 0;
	double * result = new double[6];

	///set boundary values
	result[0] = HUGE_VAL;
	result[1] = -HUGE_VAL;
	result[2] = HUGE_VAL;
	result[3] = -HUGE_VAL;
	result[4] = HUGE_VAL;
	result[5] = -HUGE_VAL;

	for(i = 0; i<size(); i++){
		Primative * temp = (Primative *) mapsto_set(objs, i);
		double * tempBounds = temp->getBounds();
		if(tempBounds[0] < result[0]){ result[0] = tempBounds[0]; }
		if(tempBounds[1] > result[1]){ result[1] = tempBounds[1]; }
		if(tempBounds[2] < result[2]){ result[2] = tempBounds[2]; }
		if(tempBounds[3] > result[3]){ result[3] = tempBounds[3]; }
		if(tempBounds[4] < result[4]){ result[4] = tempBounds[4]; }
		if(tempBounds[5] > result[5]){ result[5] = tempBounds[5]; }
		delete[](tempBounds);
	}

	return result;
}

///find the closest primative, then compute the normal
double * primativesSet::normal(double * pt)
{
	int i = 0;
	double minDistance = HUGE_VAL;
	int closestPrimative = -1;
	double tempDist;
	
	for(i = 0; i<size(); i++){
		Primative * temp = (Primative *) mapsto_set(objs, i);
		tempDist = temp->distance(pt);
		///it is assumed that we should get fairly close because of approximation.
		if(tempDist < minDistance){
			minDistance = tempDist;
			closestPrimative = i;
		}
	}
	
	Primative * temp = (Primative *) mapsto_set(objs, closestPrimative);
	double * result = temp->normal(pt);
	
	return result;
}


void primativesSet::toString()
{
	int i = 0;
	for(i = 0; i<size(); i++){
		Primative * tempPrim = getPrim(i);
		if(tempPrim->primType == PRIM_TYPE_SPHERE){
			printf("SPHERE: %f %f %f %f\n", 
				((PrimSphere*)tempPrim)->origin[0], 
				((PrimSphere*)tempPrim)->origin[1], 
				((PrimSphere*)tempPrim)->origin[2], 
				((PrimSphere*)tempPrim)->radius);
		}
		if(tempPrim->primType == PRIM_TYPE_CUBE){
			printf("CUBE: %f %f %f %f\n", 
				((PrimCube*)tempPrim)->center[0], 
				((PrimCube*)tempPrim)->center[1], 
				((PrimCube*)tempPrim)->center[2], 
				((PrimCube*)tempPrim)->size);
		}
	}
	


}


void primativesSet::constructHashTable()
{
	//construct LatticeHash
	if( hash == NULL ){
		hash = new LatticeHash();
		hash->setup(this, 2.0);
	}

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
///wrapper parser calls probeGen after generating a set of spheres from a pdb file
SurfaceObject * ProbeGenWrapper(char * ligFile, double Sphere_Radius, double resolution)
{
	FILE * file = fopen(ligFile, "r");
	if(file == NULL){ printf("Ligand File [%s] Not Found. Exitting\n", ligFile); exit(1); }

	///generate the input container.
	primativesSet * primSet = new primativesSet();

	///Set up parsing buffers.
	char * line = new char[1000];
	char * temp = new char[100];

	///get the points out of the file
	while( fgets(line, NUMCOLUMNS, file) != NULL){
		if( (line[0]=='A' && line[1]=='T' && line[2]=='O' && line[3]=='M') || 
			(line[0]=='H' && line[1]=='E' && line[2]=='T' && line[3]=='A' && line[4]=='T' && line[5]=='M') )
		{
			double * c = new double[3];
			temp[0]=line[30];	temp[1]=line[31];	temp[2]=line[32];	temp[3]=line[33];	temp[4]=line[34];	temp[5]=line[35];	temp[6]=line[36];	temp[7]=line[37];	temp[8]='\0';
			c[0] = atof(temp); 
			temp[0]=line[38];	temp[1]=line[39];	temp[2]=line[40];	temp[3]=line[41];	temp[4]=line[42];	temp[5]=line[43];	temp[6]=line[44];	temp[7]=line[45];	temp[8]='\0';
			c[1] = atof(temp); 
			temp[0]=line[46];	temp[1]=line[47];	temp[2]=line[48];	temp[3]=line[49];	temp[4]=line[50];	temp[5]=line[51];	temp[6]=line[52];	temp[7]=line[53];	temp[8]='\0';
			c[2] = atof(temp); 
			
			PrimSphere * tempSphere = new PrimSphere(c, Sphere_Radius);
			primSet->addObj( (Primative*) tempSphere);
			delete[](c);
		}
	}

	fclose(file);
	delete[](line);
	delete[](temp);
	
//	SurfaceObject * result = ProbeGen(primSet, resolution);
	SurfaceObject * result = convertPrimativesToSurface( primSet, resolution );
	

	delete(primSet);

	return result;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
///wrapper parser calls probeGen after generating a set of primatives from a PRIMS file
SurfaceObject * ProbeGenWrapper(char * primsFile, double resolution)
{
	FILE * file = fopen(primsFile, "r");
	if(file == NULL){ printf("Volumetric Primatives File [%s] Not Found.\n", primsFile); exit(1); }

	///generate the input container.
	primativesSet * primSet = new primativesSet();

	///Set up parsing buffers.
	char * line = new char[1000];
	char * ptr;
	
	///get the points out of the file
	while( fgets(line, NUMCOLUMNS, file) != NULL){
		if(line[0] == '#'){continue;}
		if(strlen(line) == 0){continue;}
		else{
			ptr = strtok(line, " \r\t\n");
			if(ptr == NULL){ continue; }
			if(strcmp(ptr, "SPHERE")==0){
				double * c = new double[3];
				double * r = new double[1];
				c[0] = atof(strtok(NULL, " \r\t\n"));
				c[1] = atof(strtok(NULL, " \r\t\n"));
				c[2] = atof(strtok(NULL, " \r\t\n"));
				r[0] = atof(strtok(NULL, " \r\t\n"));
				printf("Found Sphere: %f %f %f  Radius: %f\n", c[0], c[1], c[2], r[0]);
				
				PrimSphere * tempSphere = new PrimSphere(c, r[0]);
				primSet->addObj( (Primative *) tempSphere);
				delete[](c);
				delete[](r);
				continue;
			}
			
			if(strcmp(ptr, "CUBE")==0){
				double * c = new double[3];
				double * r = new double[1];
				c[0] = atof(strtok(NULL, " \r\t\n"));
				c[1] = atof(strtok(NULL, " \r\t\n"));
				c[2] = atof(strtok(NULL, " \r\t\n"));
				r[0] = atof(strtok(NULL, " \r\t\n"));
				printf("Found Cube: %f %f %f  halfSide: %f\n", c[0], c[1], c[2], r[0]);
				
				PrimCube * tempCube = new PrimCube(c, r[0]);
				primSet->addObj( (Primative *) tempCube);
				delete[](c);
				delete[](r);
				continue;
			}
		}
	}	
	
//	SurfaceObject * result = ProbeGen(primSet, resolution);
	SurfaceObject * result = convertPrimativesToSurface( primSet, resolution );

	delete[](line);
	delete(primSet);

	return result;
}








//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////
///This function separates and generates surface data for each nonintersecting primative,
///then uses Marching Cubes to surface the intersecting primatives.
///
///Currently this function does not procedurally generate Sphere geometry, because procedural
///generation of good spheres without retarded triangles (which will probably be a problem for my other code)
///is difficult.  However, nonintersecting cubes are procedurally generated.
///
///Note also that this does not isolate connected components of interacting primatives, which
///in the future could be separated for more rapid matching cubes operations.
SurfaceObject * convertPrimativesToSurface( primativesSet * primSet, double resolution )
{
	int i = 0;
	int j = 0;
	bool * isAlone = new bool[ primSet->size() ];
	for(i = 0; i<primSet->size(); i++){ isAlone[i] = true; }	///default to true.

	///determine which primatives intersect with each other.
	///this is a dumb pairwise test; it is assumed we dont have lots of primatives.	
	for(i = 0; i<primSet->size(); i++){
		Primative * p1 = primSet->getPrim(i);

		///preemptively ignore anything that is a sphere.
		if(p1->primType == PRIM_TYPE_SPHERE){ isAlone[i] = false; continue; }

		///test everyone that isnt yourself, and if your intersect, mark it false;
		for(j = 0; j<primSet->size(); j++){
			if(j == i){ continue; }
			else{
				Primative * p2 = primSet->getPrim(j);
				if( p1->intersects( p2 ) ){
					isAlone[i] = false;;
					break;
				}
			}	
		}
	}

	///generate a primSet 
	set_t aloneSurfaces = alloc_set(SP_MAP);
	for(i = 0; i<primSet->size(); i++){
		Primative * p = primSet->getPrim(i);
		if(isAlone[i] == true){
			SurfaceObject * tempSurf = p->genSurface();
			aloneSurfaces = associate_set(aloneSurfaces, size_set(aloneSurfaces), tempSurf);
		}
	}

	///count the number of interacting primatives.
	int counter = 0;
	for(i = 0; i<primSet->size(); i++){
//		Primative * p = primSet->getPrim(i);
		if(isAlone[i] == false){ counter++; }
	}

	primativesSet * interactingSet = NULL;
	SurfaceObject * combinedObjs = NULL;
	if(counter > 0){
		///generate a primSet of interacting primatives.
		interactingSet = new primativesSet();
		for(i = 0; i<primSet->size(); i++){
			Primative * p = primSet->getPrim(i);
			if(isAlone[i] == false){
				Primative * newPrim = p->copy();
				interactingSet->addObj(newPrim);
			}
		}
		combinedObjs = convertPrimsMarchingCubes( interactingSet, resolution );
	}	

	///combine all the individual surfaces.
	SurfaceObject * result = new SurfaceObject();
	for(i = 0; i<size_set(aloneSurfaces); i++){
		SurfaceObject * s = (SurfaceObject *) mapsto_set(aloneSurfaces, i);
		result->addObject(s);
		delete(s);			//dispose of afterwards.
	}
	free_set(aloneSurfaces);

	///add in the combined Objects
	if( combinedObjs != NULL ){
//		printf("WTF!?\n");
//		combinedObjs->toString();
		result->addObject( combinedObjs );
		delete( combinedObjs );
		delete( interactingSet );
	}

	delete[](isAlone);
	
	return result;
}












//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
////It is assumed that these prims are all co-intersecting in some manner, so we use marching
////cubes to surface them
SurfaceObject * convertPrimsMarchingCubes( primativesSet * primSet, double resolution )
{
	///get boundary values
	double xneg, xpos, yneg, ypos, zneg, zpos;
	double res = resolution;
	
	printf("Running marching Cubes to surface primatives\n");

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	///generate requisite data structures
	CubeTable * table = new CubeTable();
//	table->readFile();

	///initialize counters
	int i = 0;	int j = 0;	int k = 0;	int l = 0;  //int m = 0;

//	for(i = 0; i<size_set(coords); i++){
//		double * coord = (double *) mapsto_set(coords, i);
//		printf("%f %f %f rad: %f \n", coord[0], coord[1], coord[2], radii[i]);
//	}

	printf("PRIMSET SIZE: %i\n", primSet->size());

	double * dimensions = primSet->getPrimSetBounds();
	xneg = dimensions[0];
	xpos = dimensions[1];
	yneg = dimensions[2];
	ypos = dimensions[3];
	zneg = dimensions[4];
	zpos = dimensions[5];
	delete[](dimensions);
	
	//set up the hash table for rapid intersection testing.
	primSet->constructHashTable();

	///random nudge to offset axis aligned geometry.
	double xnudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double ynudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double znudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	xneg -= xnudge;
	yneg -= ynudge;
	zneg -= znudge;

	////////////////////////////////////////
	printf("PRIMSET DIMENSIONS:  xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	////////////////////////////////////////

	//padding on low side.
	xneg -= res;	yneg -= res;	zneg -= res;

	///+1 for fractional cube round up (int conversion truncates)
	///+1 for high side padding
	///= +2 total additions.
	int xdim = ((int) ((xpos-xneg)/res)) + 2; 	///this is the total number of CUBES.  Not lines.
	int ydim = ((int) ((ypos-yneg)/res)) + 2;
	int zdim = ((int) ((zpos-zneg)/res)) + 2;
	///Reset the high dimensions based on the sizes of the cubes, since there is fractional roundup.
	xpos = xneg + xdim*res;
	ypos = yneg + ydim*res;
	zpos = zneg + zdim*res;
	
	////////////////////////////////////////
	printf("EXPANDED DIMENSIONS: xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	printf("CUBE DIMENSIONS:     xdim: %i, ydim: %i, zdim: %i\n", xdim, ydim, zdim);
	////////////////////////////////////////

	///for each cube, we store 12 edges. (index to vertex, if there is one, -1 otherwise)
	///for each cube, we store 8 points. (inside (true) or outside (false))
	printf("Generating Data Structures...\n");
	set_t cubes = alloc_set(SP_MAP);
	set_t corners = alloc_set(SP_MAP);
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				int * theseEdges = new int[12];
				bool * thesePts = new bool[8];
				for(l = 0; l<12; l++){ theseEdges[l] = -1; }
				for(l = 0; l<8; l++){ thesePts[l] = false; }
				int currentCubeIndex = (i*ydim*zdim) + (j*zdim) + k;
				cubes = associate_set(cubes, currentCubeIndex, (ptr_t) theseEdges);
				corners = associate_set(corners, currentCubeIndex, (ptr_t) thesePts);
			}
		}
	}

	///For each cube, detect if an intersection exists on the cube edges, 
	///then store vertex indices on the cube edges.
	int progress_cube = 0;
	int progress_max = xdim*ydim*zdim;
	printf("Computing Probe Surface...\n"); fflush(stdout);
	int pointCounter = 0;
	set_t points = alloc_set(SP_MAP);
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				stdoutProgressBar(progress_cube, progress_max);

				///Allocate Stuff we will use.
//				int oldPtCounter = pointCounter;
				int currentCubeIndex = (i*ydim*zdim) + (j*zdim) + k;
				int * theseEdges = (int *) mapsto_set(cubes, currentCubeIndex);
				bool * thesePts = (bool *) mapsto_set(corners, currentCubeIndex);
				double * dist = new double[1];
				double * thisPoint = NULL;				

				///Generate the points for this cube.				
				double * c0 = new double[3]; c0[0] = xneg+(i+0)*res; c0[1] = yneg+(j+0)*res; c0[2] = zneg+(k+0)*res;
				double * c1 = new double[3]; c1[0] = xneg+(i+0)*res; c1[1] = yneg+(j+0)*res; c1[2] = zneg+(k+1)*res;
				double * c2 = new double[3]; c2[0] = xneg+(i+0)*res; c2[1] = yneg+(j+1)*res; c2[2] = zneg+(k+0)*res;
				double * c3 = new double[3]; c3[0] = xneg+(i+0)*res; c3[1] = yneg+(j+1)*res; c3[2] = zneg+(k+1)*res;
				double * c4 = new double[3]; c4[0] = xneg+(i+1)*res; c4[1] = yneg+(j+0)*res; c4[2] = zneg+(k+0)*res;
				double * c5 = new double[3]; c5[0] = xneg+(i+1)*res; c5[1] = yneg+(j+0)*res; c5[2] = zneg+(k+1)*res;
				double * c6 = new double[3]; c6[0] = xneg+(i+1)*res; c6[1] = yneg+(j+1)*res; c6[2] = zneg+(k+0)*res;
				double * c7 = new double[3]; c7[0] = xneg+(i+1)*res; c7[1] = yneg+(j+1)*res; c7[2] = zneg+(k+1)*res;

				///these bools say if youa re inside or outside
				bool d0 = primSet->isInside(c0);  
				bool d1 = primSet->isInside(c1);
				bool d2 = primSet->isInside(c2);
				bool d3 = primSet->isInside(c3);
				bool d4 = primSet->isInside(c4);
				bool d5 = primSet->isInside(c5);
				bool d6 = primSet->isInside(c6);
				bool d7 = primSet->isInside(c7);

				///Test if the points fall inside or outside the region.
				bool p0 = false, p1 = false, p2 = false, p3 = false, p4 = false, p5 = false, p6 = false, p7 = false;
				if( d0 ){ p0 = true; }   if(p0){ thesePts[0] = true;}  // printf("pt 0 is inside\n"); }          
				if( d1 ){ p1 = true; }   if(p1){ thesePts[1] = true;}  // printf("pt 1 is inside\n"); }          
				if( d2 ){ p2 = true; }   if(p2){ thesePts[2] = true;}  // printf("pt 2 is inside\n"); }          
				if( d3 ){ p3 = true; }   if(p3){ thesePts[3] = true;}  // printf("pt 3 is inside\n"); }          
				if( d4 ){ p4 = true; }   if(p4){ thesePts[4] = true;}  // printf("pt 4 is inside\n"); }          
				if( d5 ){ p5 = true; }   if(p5){ thesePts[5] = true;}  // printf("pt 5 is inside\n"); }          
				if( d6 ){ p6 = true; }   if(p6){ thesePts[6] = true;}  // printf("pt 6 is inside\n"); }          
				if( d7 ){ p7 = true; }   if(p7){ thesePts[7] = true;}  // printf("pt 7 is inside\n"); }          


				///SET THE EDGES BASED ON THE INSIDE/OUTSIDE PARITY
				///EDGE 0//////////////////////////////////////////////////////////////////////////////////////
				if((p0 && !p1) || (!p0 && p1)){
					thisPoint = primSet->intersectionPt(c0, c1);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[0] = pointCounter; pointCounter++;
				}
				///EDGE 1//////////////////////////////////////////////////////////////////////////////////////
				if((p0 && !p2) || (!p0 && p2)){
					thisPoint = primSet->intersectionPt(c0, c2);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[1] = pointCounter; pointCounter++;
				}
				///EDGE 4//////////////////////////////////////////////////////////////////////////////////////
				if((p0 && !p4) || (!p0 && p4)){
					thisPoint = primSet->intersectionPt(c0, c4);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[4] = pointCounter; pointCounter++;
				}
				////If you are on one of the high X side, get the high X edge
				if(i == xdim-1){
				///EDGE 8//////////////////////////////////////////////////////////////////////////////////////
					if((p4 && !p5) || (!p4 && p5)){
						thisPoint = primSet->intersectionPt(c4, c5);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[8] = pointCounter; pointCounter++;
					}
				///EDGE 9//////////////////////////////////////////////////////////////////////////////////////
					if((p4 && !p6) || (!p4 && p6)){
						thisPoint = primSet->intersectionPt(c4, c6);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[9] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Y side, get the high Y edge
				if(j == ydim-1){
				///EDGE 3//////////////////////////////////////////////////////////////////////////////////////
					if((p2 && !p3) || (!p2 && p3)){
						thisPoint = primSet->intersectionPt(c2, c3);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[3] = pointCounter; pointCounter++;
					}
				///EDGE 6//////////////////////////////////////////////////////////////////////////////////////
					if((p2 && !p6) || (!p2 && p6)){
						thisPoint = primSet->intersectionPt(c2, c6);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[6] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Z side, get the high Z edge
				if(k == zdim-1){
				///EDGE 2//////////////////////////////////////////////////////////////////////////////////////
					if((p1 && !p3) || (!p1 && p3)){
						thisPoint = primSet->intersectionPt(c1, c3);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[2] = pointCounter; pointCounter++;
					}
				///EDGE 5//////////////////////////////////////////////////////////////////////////////////////
					if((p1 && !p5) || (!p1 && p5)){
						thisPoint = primSet->intersectionPt(c1, c5);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[5] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high X and Y sides, get the high X and Y edge
				if(i == xdim-1 && j == ydim-1){
				///EDGE 11//////////////////////////////////////////////////////////////////////////////////////
					if((p6 && !p7) || (!p6 && p7)){
						thisPoint = primSet->intersectionPt(c6, c7);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[11] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Y and Z sides, get the high Y and Z edge
				if(j == ydim-1 && k == zdim-1){
				///EDGE 7//////////////////////////////////////////////////////////////////////////////////////
					if((p3 && !p7) || (!p3 && p7)){
						thisPoint = primSet->intersectionPt(c3, c7);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[7] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high X and Z sides, get the high X and Z edge
				if(k == zdim-1 && i == xdim-1){
				///EDGE 10//////////////////////////////////////////////////////////////////////////////////////
					if((p5 && !p7) || (!p5 && p7)){
						thisPoint = primSet->intersectionPt(c5, c7);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[10] = pointCounter; pointCounter++;
					}
				}

//				///RESET EARLIER EDGES BECAUSE WE DID NOT SET THEM BEFORE
				int * tempCube;
				if(theseEdges[0] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[8] = theseEdges[0];
					}
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[3] = theseEdges[0];
					}
					if(i != 0 && j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[11] = theseEdges[0];
					}
				}
				if(theseEdges[1] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[9] = theseEdges[1];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[2] = theseEdges[1];
					}
					if(i != 0 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[10] = theseEdges[1];
					}
				}
				if(theseEdges[2] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[10] = theseEdges[2];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[1] = theseEdges[2];
					}
					if(i != 0 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[9] = theseEdges[2];
					}
				}
				if(theseEdges[3] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[11] = theseEdges[3];
					}
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[0] = theseEdges[3];
					}
					if(i != 0 && j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[8] = theseEdges[3];
					}
				}
				if(theseEdges[4] != -1){
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[6] = theseEdges[4];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[5] = theseEdges[4];
					}
					if(j != 0 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k-1);
						tempCube[7] = theseEdges[4];
					}
				}
				if(theseEdges[5] != -1){
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[7] = theseEdges[5];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[4] = theseEdges[5];
					}
					if(j != 0 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+1);
						tempCube[6] = theseEdges[5];
					}
				}
				if(theseEdges[6] != -1){
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[4] = theseEdges[6];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[7] = theseEdges[6];
					}
					if(j != ydim-1 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k-1);
						tempCube[5] = theseEdges[6];
					}
				}
				if(theseEdges[7] != -1){
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[5] = theseEdges[7];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[6] = theseEdges[7];
					}
					if(j != ydim-1 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+1);
						tempCube[4] = theseEdges[7];
					}
				}
				if(theseEdges[8] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[0] = theseEdges[8];
					}
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[11] = theseEdges[8];
					}
					if(i != xdim-1 && j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[3] = theseEdges[8];
					}
				}
				if(theseEdges[9] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[1] = theseEdges[9];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[10] = theseEdges[9];
					}
					if(i != xdim-1 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[2] = theseEdges[9];
					}
				}
				if(theseEdges[10] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[2] = theseEdges[10];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[9] = theseEdges[10];
					}
					if(i != xdim-1 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[1] = theseEdges[10];
					}
				}
				if(theseEdges[11] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[3] = theseEdges[11];
					}
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[8] = theseEdges[11];
					}
					if(i != xdim-1 && j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[0] = theseEdges[11];
					}
				}

				delete[](dist);
				delete[](c0); delete[](c1); delete[](c2); delete[](c3);
				delete[](c4); delete[](c5); delete[](c6); delete[](c7);
				progress_cube++;
			}
		}		
	}

	///printf("POINTS GENERATED: %i\n", size_set(points));

	////Second Pass: generate Triangle List.
	printf("Generating triangles...\n"); fflush(stdout);
	int triListSize = 0;
	int triListCap = 100;
	int ** triList = new int*[triListCap];
	double ** normList = new double*[triListCap];
	vertTriSet * triSet = new vertTriSet();

	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				///look up the internal/external points, as stored before.
				int currentCubeIndex = (i*ydim*zdim) + (j*zdim) + k;
				bool * thesePts = (bool *) mapsto_set(corners, currentCubeIndex);
				int p0i = 0; if(thesePts[0]){ p0i = 1; }
				int p1i = 0; if(thesePts[1]){ p1i = 1; }
				int p2i = 0; if(thesePts[2]){ p2i = 1; }
				int p3i = 0; if(thesePts[3]){ p3i = 1; }
				int p4i = 0; if(thesePts[4]){ p4i = 1; }
				int p5i = 0; if(thesePts[5]){ p5i = 1; }
				int p6i = 0; if(thesePts[6]){ p6i = 1; }
				int p7i = 0; if(thesePts[7]){ p7i = 1; }

				///figure out what case we in the cubetable we will use.
				int cubeTableIndex = p0i + 2*p1i + 4*p2i + 8*p3i + 16*p4i + 32*p5i + 64*p6i + 128*p7i;
				int numPolysToAdd = table->sizes[cubeTableIndex];

				///create triangles for this cube.
				if(numPolysToAdd > 0){
					////expand triList and normList if they need to be expaded.
					if(triListSize + numPolysToAdd > triListCap){
						int ** newTriList = new int*[triListCap*2];
						double ** newNormList = new double*[triListCap*2];
						for(l = 0; l<triListSize; l++){
							newTriList[l] = triList[l];
							newNormList[l] = normList[l];
						}
						delete[](triList);
						delete[](normList);
						triList = newTriList;
						normList = newNormList;
						triListCap = triListCap*2;
					}
					///now create the triangles
//					printf("Creating %i Triangles.  CubeTableIndex: %i\n", numPolysToAdd, cubeTableIndex);
					int * polyList = table->polys[cubeTableIndex];
					int * thisCube = (int *) mapsto_set(cubes, currentCubeIndex);
					for(l = 0; l<numPolysToAdd; l++){
						int * newTriangle = new int[3];

						//printf("getting the point indices: numPolys: %i Edges: [%i %i %i]\n", numPolysToAdd, polyList[0], polyList[1], polyList[2] );
						///add the triangle into the triList.
						triList[triListSize] = newTriangle;

						///reversed order here to reverse norms - original polyList is reversed.
						newTriangle[0] = thisCube[polyList[3*l+0]];
						newTriangle[1] = thisCube[polyList[3*l+2]];	
						newTriangle[2] = thisCube[polyList[3*l+1]];

						//printf("adding indices to the triList\n");
						//store point indices into triList, for this triangle
						triSet->addEdge(newTriangle[0], triListSize);	
						triSet->addEdge(newTriangle[1], triListSize);
						triSet->addEdge(newTriangle[2], triListSize);

						//printf("Getting point coords: newTriangle[0]: %i, newTriangle[1]: %i newTriangle[2]: %i\n", newTriangle[0], newTriangle[1], newTriangle[2]);
						//get point Vals (do not delete these)
						double * firstPt = (double *) mapsto_set(points, newTriangle[0]);
						double * secndPt = (double *) mapsto_set(points, newTriangle[1]);
						double * thirdPt = (double *) mapsto_set(points, newTriangle[2]);

						///now get the triangle normal
						double * normal;
						
						///if the area is nonzero, proceed as normal
						if(  vectorSize(firstPt[0]-secndPt[0], firstPt[1]-secndPt[1], firstPt[2]-secndPt[2]) > 0 &&
							vectorSize(secndPt[0]-thirdPt[0], secndPt[1]-thirdPt[1], secndPt[2]-thirdPt[2]) > 0 &&
							vectorSize(firstPt[0]-thirdPt[0], firstPt[1]-thirdPt[1], firstPt[2]-thirdPt[2]) > 0 )
						{
							///get normals (delete first, second)
							//printf("Calcing first diff for X-prod\n");
							double * first = new double[3];
							first[0] = (secndPt[0] - firstPt[0]);
							first[1] = (secndPt[1] - firstPt[1]);
							first[2] = (secndPt[2] - firstPt[2]);
	
							//printf("Calcing 2nd diff for X-prod\n");
							double * second = new double[3];
							second[0] = (thirdPt[0] - firstPt[0]);
							second[1] = (thirdPt[1] - firstPt[1]);
							second[2] = (thirdPt[2] - firstPt[2]);
	
							//printf("Calcing X-prod\n");
							normal = crossProd(first, second);
							delete[](first);
							delete[](second);
						}
						//if the area is zero, get the normal based on the primative centers.
						else{
							double * avg = new double[3];
							avg[0] = (firstPt[0] + secndPt[0] + thirdPt[0])/3;
							avg[1] = (firstPt[1] + secndPt[1] + thirdPt[1])/3;
							avg[2] = (firstPt[2] + secndPt[2] + thirdPt[2])/3;

							normal = primSet->normal(avg);

							delete[](avg);
						}
						double * normalizedNormal = normalizeVector(normal);
						normList[triListSize] = normalizedNormal;

						////clear and increment
						delete[](normal);

						triListSize++;
					}
				}
				else{
					///commented out bc it pritns out too much
					//printf("No Tris To Make!\n");
				}
				
			}
		}
	}

	////Third Pass: produce actual triangles and normals.
	printf("Generating Triangles and Normals\n"); fflush(stdout);
	int numPts;
	double * pts;
	double * norms;
	int numTris;
	int * tris;
	double * triangleNorms;

	double reverseNorm = 1;

	///Copy over the points	
	numPts = size_set(points);
	pts = new double[3*size_set(points)];
	for(i = 0; i<size_set(points); i++){
		double * tempPt = (double *) mapsto_set(points, i);
		pts[3*i+0] = tempPt[0];	pts[3*i+1] = tempPt[1];	pts[3*i+2] = tempPt[2];
	}
	///Copy over the triangles and triangle normals
	numTris = triListSize;
	tris = new int[3*triListSize];
	triangleNorms = new double[3*triListSize];
	for(i = 0; i<triListSize; i++){
		tris[3*i+0] = triList[i][0];	
		tris[3*i+1] = triList[i][1];	
		tris[3*i+2] = triList[i][2];
		triangleNorms[3*i+0] = reverseNorm*normList[i][0];	
		triangleNorms[3*i+1] = reverseNorm*normList[i][1];	
		triangleNorms[3*i+2] = reverseNorm*normList[i][2];
	}
	
	////STEP 3: Generate Averaged normals
	norms = new double[3*numPts];///this is output
	for(i = 0; i<triListSize; i++){
		int * thisTri = triList[i];
		////for each vertex in the triangle, get the list of adjacent triangles, and get the normals.
		///This loop processes each point in the triangle separately.
		for(j = 0; j<3; j++){
			int thisPt = thisTri[j];
//			double * thisVector = (double *) mapsto_set(points, thisPt);  ///this gets the coordinates of the point
			set_t neighbors = triSet->getNeighbors(thisPt); ///this gets the triangles that touch this point.
			int numNeigh = size_set(neighbors);  ///this is the number of triangles adjacent
			double * normal = new double[3];	///the normal for this point
			normal[0] = 0;		normal[1] = 0;		normal[2] = 0;
			for(k = 0; k<numNeigh; k++){
				normal[0] = normal[0] + normList[neighbors[k]][0];	///add up the normals of the adjacent traignels
				normal[1] = normal[1] + normList[neighbors[k]][1];
				normal[2] = normal[2] + normList[neighbors[k]][2];
//				if( 	normList[neighbors[k]][0] != normList[neighbors[k]][0] || 
//					normList[neighbors[k]][1] != normList[neighbors[k]][1] || 
//					normList[neighbors[k]][2] != normList[neighbors[k]][2] ){ printf("OMFG!\n"); }
			}
			normal[0] /= (double) numNeigh;	///average the normals of adjacent triangles
			normal[1] /= (double) numNeigh;
			normal[2] /= (double) numNeigh;
//			if( normal[0] != normal[0] || normal[1] != normal[1] || normal[2] != normal[2] ){ printf("OMFG!! numNeigh: %f\n", (double) numNeigh); }
			
			double * normedNormal = normalizeVector(normal);   //This is now the averaged normal; normalize it.
			norms[3*thisPt+0] = reverseNorm*normedNormal[0];	///Store the result into norms for this point in the triangle.
			norms[3*thisPt+1] = reverseNorm*normedNormal[1];
			norms[3*thisPt+2] = reverseNorm*normedNormal[2];
//			if( norms[3*thisPt+0]!=norms[3*thisPt+0] || norms[3*thisPt+1]!=norms[3*thisPt+1] || norms[3*thisPt+2]!=norms[3*thisPt+2] ){ printf("OMFG!!!\n"); }

			free_set(neighbors);
			delete[](normedNormal);
			delete[](normal);
		}
	}



	////STEP 4: Clean up all the memory allocations
	///////////////////////////////////////////////////////////////////////////////
	for(i = 0; i<size_set(points); i++){ double * vector = (double *) mapsto_set(points, i); delete[](vector); }
	for(i = 0; i<size_set(corners); i++){ int * corner = (int *) mapsto_set(corners, i); delete[](corner); }
	for(i = 0; i<size_set(cubes); i++){ int * cube = (int *) mapsto_set(cubes, i); delete[](cube); }
	for(i = 0; i<triListSize; i++){ delete[](triList[i]); delete[](normList[i]); }
	free_set(points);
	free_set(corners);
	free_set(cubes);
	delete[](triList);
	delete[](normList);
	delete(triSet);
	///////////////////////////////////////////////////////////////////////////////
	////STEP 4: Clean up all the memory allocations

	////STEP 5: Generate the output:	
	///////////////////////////////////////////////////////////////////////////////
	SurfaceObject * results = new SurfaceObject (numPts, pts, norms, numTris, tris, triangleNorms);

	///debug
	//generateSURF(results, NULL, "TMESHTEST.SURF");
	///////////////////////////////////////////////////////////////////////////////
	////STEP 5: Generate the output:	
	

	////Step 6: Clear the output data:
	///////////////////////////////////////////////////////////////////////////////
	delete[](pts);
	delete[](norms);
	delete[](tris);
	delete[](triangleNorms);
	///////////////////////////////////////////////////////////////////////////////
	////Step 6: Clear the output data:

	delete(table);

	/////Return the result
	///////////////////////////////////////////////////////////////////////////////
	return results;



}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////
////This method replaces getAverageIntersectionPt, used in generateSimplifiedSurface
////in order to deal with more complex boolean processing.
////
///This function is assuming that it only gets called when there is definitely an intersection.
////
////  Because we are averaging the underlying geometry, we implement a policy:
////    When we can tell which answer is correct, we provide it.  otherwise we average.
double * getBooleanPrimGeometry(bool p1a, bool p1b, bool p2a, bool p2b, 
						TriangleLatticeHash * hasha, primativesSet * primSet, 
						double * p1, double * p2, int op)
{
	double * result = NULL;

	bool flipped = false;	
	if( !(p1b==p2b) ){
		//if they are not the same, and p1b is outside, we switch them to make sure 
		//that intersectionPt() always gets the first one "inside"
		if(!p1b){
			flipped = true;
		}
	}

	///first, detect if there is only one intersection.
	///This function is assuming that it only gets called when there is definitely an intersection.
	if( (p1a==p2a) && !(p1b==p2b) ){ 
		if(flipped){ result = primSet->intersectionPt(p2, p1); return result; }
		else{ result = primSet->intersectionPt(p1, p2); return result; }
	}
	if( !(p1a==p2a) && (p1b==p2b) ){ result = hasha->getAverageIntersectionPt(p1, p2); return result; }

	///this means that there are two intersections.  This is the difficult case.
	/////first, get the intersections.  The answer is one of these.
	double * pointa = hasha->getAverageIntersectionPt(p1, p2);
	double * pointb;
	if(flipped){ pointb = primSet->intersectionPt(p2, p1); }
	else{ pointb = primSet->intersectionPt(p1, p2); }
		
	double * average = new double[3];
	average[0] = (pointa[0]+pointb[0])/2;
	average[1] = (pointa[1]+pointb[1])/2;
	average[2] = (pointa[2]+pointb[2])/2;

	//double * dist = new double[1];
	bool aInsideb = primSet->isInside(pointa);
	bool bInsidea = hasha->isInside(pointb);

	///handle each case:
	if(op == SURFACE_BOOLEAN_UNION ){
		if(aInsideb == bInsidea){ result = average; }
		if(aInsideb && !bInsidea){ result = pointb; }
		if(!aInsideb && bInsidea){ result = pointa; }
	}
	if(op == SURFACE_BOOLEAN_INTERSECT ){
		if(aInsideb == bInsidea){ result = average; }
		if(aInsideb && !bInsidea){ result = pointa; }
		if(!aInsideb && bInsidea){ result = pointb; }
	}
	if(op == SURFACE_BOOLEAN_DIFFERENCE ){
		if( aInsideb &&  bInsidea){ result = pointb; }
		if(!aInsideb && !bInsidea){ result = pointa; }
		if( aInsideb && !bInsidea){ result = pointa; }
		if(!aInsideb &&  bInsidea){ result = pointb; }
	}

	double * output = new double[3];
	output[0] = result[0];
	output[1] = result[1];
	output[2] = result[2];

	delete[](pointa);
	delete[](pointb);
	delete[](average);
	//delete[](dist);

	return output;					

}





//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
///Boolean marching cubes using a probeset.  Leverages the accuracy of the 
///probeset's analytic computations.
SurfaceObject * booleanProbeProcessing(SurfaceObject * surf1, bool neg1, primativesSet * primSet, bool neg2, int op, double res)
{
	printf("Running Analytical Boolean Probe Operation.");
	printf("Note that this has NOT been tested on inside out probes.");

	///generate requisite data structures
	double hashRes = 2.0;

	SurfaceObject * result;
	///detect the degenerate case - both are empty.
	if( (surf1->numPoints==0) && (primSet->size()==0) ){
		result = new SurfaceObject();
		return result;
	}

	///detect the semi-degenerate case.
	if( (surf1->numPoints==0) || (primSet->size()==0) ){
		double rep1 = 1;
		double rep2 = 2;
		double rep3;
		
		if(op == SURFACE_BOOLEAN_DIFFERENCE){
			neg2 = !neg2;
			op = SURFACE_BOOLEAN_INTERSECT;
		}	
		
		if(surf1->numPoints==0){ rep1 = 0; }
		if(primSet->size()==0){ rep2 = 0; }
			
		if(neg1 && surf1->numPoints==0){ rep1 = HUGE_VAL; }
		if(neg2 && primSet->size()==0){ rep2 = HUGE_VAL; }
			
		if(neg1 && surf1->numPoints!=0){ rep1 = -1; }
		if(neg2 && primSet->size()!=0){ rep2 = -2; }
		
		if(op == SURFACE_BOOLEAN_UNION){
			rep3 = rep1+rep2;
		}
		if(op == SURFACE_BOOLEAN_INTERSECT){
			//rep3 = rep1*rep2;
			rep3 = 0;
			if(fabs(rep1) > 5){ rep3 = rep2; }
			if(fabs(rep2) > 5){ rep3 = rep1; }
			
		}

		//evaluate rep3;
		if(rep3 == 0){
			result = new SurfaceObject();
		}
		if(rep3 == -1){
			result = surf1->copy(NULL); result->flipNormals();
		}
		if(rep3 == -2){
			result = convertPrimativesToSurface( primSet, res );
			result->flipNormals();
		}
		if(rep3 == 1){
			result = surf1->copy(NULL);
		}
		if(rep3 == 2){
			result = convertPrimativesToSurface( primSet, res );
		}
		//note that the empty space and the full space are identical and indistinguishable.
		if( fabs(rep3) > 5 ){
			result = new SurfaceObject(); 
		}

		return result;
	}

	///initialize counters
	int i = 0;	int j = 0;	int k = 0;	int l = 0;

	///get boundary values
	double xneg = HUGE_VAL;		double xpos = -HUGE_VAL;
	double yneg = HUGE_VAL;		double ypos = -HUGE_VAL;
	double zneg = HUGE_VAL;		double zpos = -HUGE_VAL;

	///temp boundary values
	double xneg1 = HUGE_VAL;		double xpos1 = -HUGE_VAL;
	double yneg1 = HUGE_VAL;		double ypos1 = -HUGE_VAL;
	double zneg1 = HUGE_VAL;		double zpos1 = -HUGE_VAL;
	///temp boundary values
	double * tempPrimSetBounds = primSet->getPrimSetBounds();
	double xneg2 = tempPrimSetBounds[0];		double xpos2 = tempPrimSetBounds[1];
	double yneg2 = tempPrimSetBounds[2];		double ypos2 = tempPrimSetBounds[3];
	double zneg2 = tempPrimSetBounds[4];		double zpos2 = tempPrimSetBounds[5];
	delete[](tempPrimSetBounds);

	for(i = 0; i<surf1->numPoints; i++){
		if(surf1->surfacePoints[3*i+0] < xneg1){ xneg1 = surf1->surfacePoints[3*i+0]; }
		if(surf1->surfacePoints[3*i+0] > xpos1){ xpos1 = surf1->surfacePoints[3*i+0]; }
		if(surf1->surfacePoints[3*i+1] < yneg1){ yneg1 = surf1->surfacePoints[3*i+1]; }
		if(surf1->surfacePoints[3*i+1] > ypos1){ ypos1 = surf1->surfacePoints[3*i+1]; }
		if(surf1->surfacePoints[3*i+2] < zneg1){ zneg1 = surf1->surfacePoints[3*i+2]; }
		if(surf1->surfacePoints[3*i+2] > zpos1){ zpos1 = surf1->surfacePoints[3*i+2]; }
	}

	////////////////////////////////////////
	printf("SURF1 DIMENSIONS:  xneg1: %f xpos1: %f yneg1: %f ypos1: %f zneg1: %f zpos1: %f\n", xneg1, xpos1, yneg1, ypos1, zneg1, zpos1);
	printf("PRIMSET DIMENSIONS:  xneg2: %f xpos2: %f yneg2: %f ypos2: %f zneg2: %f zpos2: %f\n", xneg2, xpos2, yneg2, ypos2, zneg2, zpos2);
	////////////////////////////////////////

	if(op == SURFACE_BOOLEAN_UNION ){
		if(xneg1 < xneg2){ xneg = xneg1; }	else{ xneg = xneg2; }
		if(xpos1 > xpos2){ xpos = xpos1; }	else{ xpos = xpos2; }
		if(yneg1 < yneg2){ yneg = yneg1; }	else{ yneg = yneg2; }
		if(ypos1 > ypos2){ ypos = ypos1; }	else{ ypos = ypos2; }
		if(zneg1 < zneg2){ zneg = zneg1; }	else{ zneg = zneg2; }
		if(zpos1 > zpos2){ zpos = zpos1; }	else{ zpos = zpos2; }
	}
	if(op == SURFACE_BOOLEAN_INTERSECT ){
		if(xneg1 < xneg2){ xneg = xneg2; }	else{ xneg = xneg1; }
		if(xpos1 > xpos2){ xpos = xpos2; }	else{ xpos = xpos1; }
		if(yneg1 < yneg2){ yneg = yneg2; }	else{ yneg = yneg1; }
		if(ypos1 > ypos2){ ypos = ypos2; }	else{ ypos = ypos1; }
		if(zneg1 < zneg2){ zneg = zneg2; }	else{ zneg = zneg1; }
		if(zpos1 > zpos2){ zpos = zpos2; }	else{ zpos = zpos1; }
	}
	
	if(op == SURFACE_BOOLEAN_DIFFERENCE ){
		xneg = xneg1;
		xpos = xpos1;
		yneg = yneg1;
		ypos = ypos1;
		zneg = zneg1;
		zpos = zpos1;
	}

	///if the surface provided is trivial/nonexistant, return the non-surface for processing
	///This will reduce problems later, because we handle non-surfaces correctly.
	if( (xneg>xpos) || (yneg>ypos) || (zneg>zpos) ){
		SurfaceObject * result = new SurfaceObject();
		return result;
	}

	//padding on low side.
	xneg -= res;	yneg -= res;	zneg -= res;

	xneg = floor(xneg);
	yneg = floor(yneg);
	zneg = floor(zneg);
	xpos = ceil(xpos);
	ypos = ceil(ypos);
	zpos = ceil(zpos);

	///random nudge to offset axis aligned geometry.
	double xnudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double ynudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double znudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	xneg -= xnudge;
	yneg -= ynudge;
	zneg -= znudge;

	////////////////////////////////////////
	printf("SURFACE DIMENSIONS:  xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos); fflush(stdout);
	////////////////////////////////////////

	///+1 for fractional cube round up (int conversion truncates)
	///+1 for high side padding
	///= +2 total additions.
	int xdim = ((int) ((xpos-xneg)/res)) + 2; 	///this is the total number of CUBES.  Not lines.
	int ydim = ((int) ((ypos-yneg)/res)) + 2;
	int zdim = ((int) ((zpos-zneg)/res)) + 2;
	///Reset the high dimensions based on the sizes of the cubes, since there is fractional roundup.
	xpos = xneg + xdim*res;
	ypos = yneg + ydim*res;
	zpos = zneg + zdim*res;
	
	////////////////////////////////////////
	printf("EXPANDED DIMENSIONS: xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos); fflush(stdout);
	printf("CUBE DIMENSIONS:     xdim: %i, ydim: %i, zdim: %i\n", xdim, ydim, zdim); fflush(stdout);
	////////////////////////////////////////

	///set up the hashing now that the dimensions are checked, and we have not exitted.
	TriangleLatticeHash * hash1 = new TriangleLatticeHash(surf1, hashRes, "BooleanSurface1");
	hash1->isInsideOut = neg1;

//we dont need this for primset computations
//	TriangleLatticeHash * hash2 = new TriangleLatticeHash(surf2, hashRes, "BooleanSurface2");		
	primSet->isInsideOut = neg2;	///this is true if it is negative (insideout) otherwise its false;

	if(hash1->isInsideOut){ printf("HASH IS INSIDE OUT\n"); }
	if(primSet->isInsideOut){ printf("PRIMSET IS INSIDE OUT\n"); }

	CubeTable * table = new CubeTable();
//	table->readFile();

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
#ifdef FORCE_LOGICAL_CONSISTANCY							///
	///*2 is min mathematically possible, *3 is overlap space.  ///
	///*2 prolly works, but we dont need the ram.			///
	int maxLookback = (ydim+1)*(zdim+1)*3;					///
	hash1->setCacheSize(maxLookback);						///
//	hash2->setCacheSize(maxLookback);						///
#endif												///
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

	///for each cube, we store 12 edges. (index to vertex, if there is one, -1 otherwise)
	///for each cube, we store 8 points. (inside (true) or outside (false))
	printf("Generating Data Structures...\n"); fflush(stdout);
	set_t cubes = alloc_set(SP_MAP);
	set_t corners = alloc_set(SP_MAP);
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				int * theseEdges = new int[12];
				bool * thesePts = new bool[8];
				for(l = 0; l<12; l++){ theseEdges[l] = -1; }
				for(l = 0; l<8; l++){ thesePts[l] = false; }
				int currentCubeIndex = (i*ydim*zdim) + (j*zdim) + k;
				cubes = associate_set(cubes, currentCubeIndex, (ptr_t) theseEdges);
				corners = associate_set(corners, currentCubeIndex, (ptr_t) thesePts);
			}
		}
	}

	///For each cube, detect if an intersection exists on the cube edges, 
	///then store vertex indices on the cube edges.
	int progress_cube = 0;
	int progress_max = xdim*ydim*zdim;
	printf("Computing Binary Geometry...\n"); fflush(stdout);
	int pointCounter = 0;
	set_t points = alloc_set(SP_MAP);
	for(i = 0; i<xdim; i++){
		//printf("i: %i\n", i);
		for(j = 0; j<ydim; j++){
			//printf("i: %i  j: %i\n", i, j);
			for(k = 0; k<zdim; k++){
				stdoutProgressBar(progress_cube, progress_max);

//				i = 6;	j = 18;	k = 50;
//				printf("Starting new cube!: i: %i xdim: %i j: %i ydim: %i k: %i zdim: %i\n", i, xdim, j, ydim, k, zdim);
//				if(i == 16 && j==22 && k == 23){
//					printf("WTF!\n");
//				}

				///Allocate Stuff we will use.
//				int oldPtCounter = pointCounter;
				int currentCubeIndex = (i*ydim*zdim) + (j*zdim) + k;
				int * theseEdges = (int *) mapsto_set(cubes, currentCubeIndex);
				bool * thesePts = (bool *) mapsto_set(corners, currentCubeIndex);
				//double * dist = new double[1];
				double * thisPoint = NULL;				

				///Generate the points for this cube.				
				double * c0 = new double[3]; c0[0] = xneg+(i+0)*res; c0[1] = yneg+(j+0)*res; c0[2] = zneg+(k+0)*res;
				double * c1 = new double[3]; c1[0] = xneg+(i+0)*res; c1[1] = yneg+(j+0)*res; c1[2] = zneg+(k+1)*res;
				double * c2 = new double[3]; c2[0] = xneg+(i+0)*res; c2[1] = yneg+(j+1)*res; c2[2] = zneg+(k+0)*res;
				double * c3 = new double[3]; c3[0] = xneg+(i+0)*res; c3[1] = yneg+(j+1)*res; c3[2] = zneg+(k+1)*res;
				double * c4 = new double[3]; c4[0] = xneg+(i+1)*res; c4[1] = yneg+(j+0)*res; c4[2] = zneg+(k+0)*res;
				double * c5 = new double[3]; c5[0] = xneg+(i+1)*res; c5[1] = yneg+(j+0)*res; c5[2] = zneg+(k+1)*res;
				double * c6 = new double[3]; c6[0] = xneg+(i+1)*res; c6[1] = yneg+(j+1)*res; c6[2] = zneg+(k+0)*res;
				double * c7 = new double[3]; c7[0] = xneg+(i+1)*res; c7[1] = yneg+(j+1)*res; c7[2] = zneg+(k+1)*res;

				///Test if the points fall inside or outside the region.
				bool p0a = false, p1a = false, p2a = false, p3a = false, p4a = false, p5a = false, p6a = false, p7a = false;
				int h0a = hash1->insideHighlightedCube(c0); if(h0a==HIGHLIGHT_INTERIOR){p0a=true;} if(h0a==HIGHLIGHT_EXTERIOR){p0a=false;}
				int h1a = hash1->insideHighlightedCube(c1); if(h1a==HIGHLIGHT_INTERIOR){p1a=true;} if(h1a==HIGHLIGHT_EXTERIOR){p1a=false;} 
				int h2a = hash1->insideHighlightedCube(c2); if(h2a==HIGHLIGHT_INTERIOR){p2a=true;} if(h2a==HIGHLIGHT_EXTERIOR){p2a=false;} 
				int h3a = hash1->insideHighlightedCube(c3); if(h3a==HIGHLIGHT_INTERIOR){p3a=true;} if(h3a==HIGHLIGHT_EXTERIOR){p3a=false;} 
				int h4a = hash1->insideHighlightedCube(c4); if(h4a==HIGHLIGHT_INTERIOR){p4a=true;} if(h4a==HIGHLIGHT_EXTERIOR){p4a=false;} 
				int h5a = hash1->insideHighlightedCube(c5); if(h5a==HIGHLIGHT_INTERIOR){p5a=true;} if(h5a==HIGHLIGHT_EXTERIOR){p5a=false;} 
				int h6a = hash1->insideHighlightedCube(c6); if(h6a==HIGHLIGHT_INTERIOR){p6a=true;} if(h6a==HIGHLIGHT_EXTERIOR){p6a=false;} 
				int h7a = hash1->insideHighlightedCube(c7); if(h7a==HIGHLIGHT_INTERIOR){p7a=true;} if(h7a==HIGHLIGHT_EXTERIOR){p7a=false;} 
				bool p0b = false, p1b = false, p2b = false, p3b = false, p4b = false, p5b = false, p6b = false, p7b = false;
				if( primSet->isInside(c0) ){ p0b=true; }else{ p0b=false; }
				if( primSet->isInside(c1) ){ p1b=true; }else{ p1b=false; }
				if( primSet->isInside(c2) ){ p2b=true; }else{ p2b=false; }
				if( primSet->isInside(c3) ){ p3b=true; }else{ p3b=false; }
				if( primSet->isInside(c4) ){ p4b=true; }else{ p4b=false; }
				if( primSet->isInside(c5) ){ p5b=true; }else{ p5b=false; }
				if( primSet->isInside(c6) ){ p6b=true; }else{ p6b=false; }
				if( primSet->isInside(c7) ){ p7b=true; }else{ p7b=false; }

				//use raytesting if not sure.
				if(h0a==NOT_HIGHLIGHTED){p0a = hash1->isInside( c0 ); if(p0a){ hash1->floodFillStart(c0, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c0, HIGHLIGHT_EXTERIOR);} }
				if(h1a==NOT_HIGHLIGHTED){p1a = hash1->isInside( c1 ); if(p1a){ hash1->floodFillStart(c1, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c1, HIGHLIGHT_EXTERIOR);} }
				if(h2a==NOT_HIGHLIGHTED){p2a = hash1->isInside( c2 ); if(p2a){ hash1->floodFillStart(c2, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c2, HIGHLIGHT_EXTERIOR);} }
				if(h3a==NOT_HIGHLIGHTED){p3a = hash1->isInside( c3 ); if(p3a){ hash1->floodFillStart(c3, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c3, HIGHLIGHT_EXTERIOR);} }
				if(h4a==NOT_HIGHLIGHTED){p4a = hash1->isInside( c4 ); if(p4a){ hash1->floodFillStart(c4, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c4, HIGHLIGHT_EXTERIOR);} }
				if(h5a==NOT_HIGHLIGHTED){p5a = hash1->isInside( c5 ); if(p5a){ hash1->floodFillStart(c5, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c5, HIGHLIGHT_EXTERIOR);} }
				if(h6a==NOT_HIGHLIGHTED){p6a = hash1->isInside( c6 ); if(p6a){ hash1->floodFillStart(c6, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c6, HIGHLIGHT_EXTERIOR);} }
				if(h7a==NOT_HIGHLIGHTED){p7a = hash1->isInside( c7 ); if(p7a){ hash1->floodFillStart(c7, HIGHLIGHT_INTERIOR); } else{hash1->floodFillStart(c7, HIGHLIGHT_EXTERIOR);} }

				//implement boolean arithmetic here:
				if(op == SURFACE_BOOLEAN_UNION ){
					if(p0a || p0b){ thesePts[0] = true;}
					if(p1a || p1b){ thesePts[1] = true;}
					if(p2a || p2b){ thesePts[2] = true;}
					if(p3a || p3b){ thesePts[3] = true;}
					if(p4a || p4b){ thesePts[4] = true;}
					if(p5a || p5b){ thesePts[5] = true;}
					if(p6a || p6b){ thesePts[6] = true;}
					if(p7a || p7b){ thesePts[7] = true;}
				}
				if(op == SURFACE_BOOLEAN_INTERSECT ){
					if(p0a && p0b){ thesePts[0] = true;}
					if(p1a && p1b){ thesePts[1] = true;}
					if(p2a && p2b){ thesePts[2] = true;}
					if(p3a && p3b){ thesePts[3] = true;}
					if(p4a && p4b){ thesePts[4] = true;}
					if(p5a && p5b){ thesePts[5] = true;}
					if(p6a && p6b){ thesePts[6] = true;}
					if(p7a && p7b){ thesePts[7] = true;}
				}
				if(op == SURFACE_BOOLEAN_DIFFERENCE ){
					if(p0a && !p0b){ thesePts[0] = true;}
					if(p1a && !p1b){ thesePts[1] = true;}
					if(p2a && !p2b){ thesePts[2] = true;}
					if(p3a && !p3b){ thesePts[3] = true;}
					if(p4a && !p4b){ thesePts[4] = true;}
					if(p5a && !p5b){ thesePts[5] = true;}
					if(p6a && !p6b){ thesePts[6] = true;}
					if(p7a && !p7b){ thesePts[7] = true;}
				}


				///SET THE EDGES BASED ON THE INSIDE/OUTSIDE PARITY
				///EDGE 0//////////////////////////////////////////////////////////////////////////////////////
				if((thesePts[0] && !thesePts[1]) || (!thesePts[0] && thesePts[1])){
					thisPoint = getBooleanPrimGeometry(p0a, p0b, p1a, p1b, hash1, primSet, c0, c1, op);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[0] = pointCounter; pointCounter++;
				}
				///EDGE 1//////////////////////////////////////////////////////////////////////////////////////
				if((thesePts[0] && !thesePts[2]) || (!thesePts[0] && thesePts[2])){
					thisPoint = getBooleanPrimGeometry(p0a, p0b, p2a, p2b, hash1, primSet, c0, c2, op);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[1] = pointCounter; pointCounter++;
				}
				///EDGE 4//////////////////////////////////////////////////////////////////////////////////////
				if((thesePts[0] && !thesePts[4]) || (!thesePts[0] && thesePts[4])){
					thisPoint = getBooleanPrimGeometry(p0a, p0b, p4a, p4b, hash1, primSet, c0, c4, op);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[4] = pointCounter; pointCounter++;
				}
				////If you are on one of the high X side, get the high X edge
				if(i == xdim-1){
				///EDGE 8//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[4] && !thesePts[5]) || (!thesePts[4] && thesePts[5])){
						thisPoint = getBooleanPrimGeometry(p4a, p4b, p5a, p5b, hash1, primSet, c4, c5, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[8] = pointCounter; pointCounter++;
					}
				///EDGE 9//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[4] && !thesePts[6]) || (!thesePts[4] && thesePts[6])){
						thisPoint = getBooleanPrimGeometry(p4a, p4b, p6a, p6b, hash1, primSet, c4, c6, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[9] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Y side, get the high Y edge
				if(j == ydim-1){
				///EDGE 3//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[2] && !thesePts[3]) || (!thesePts[2] && thesePts[3])){
						thisPoint = getBooleanPrimGeometry(p2a, p2b, p3a, p3b, hash1, primSet, c2, c3, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[3] = pointCounter; pointCounter++;
					}
				///EDGE 6//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[2] && !thesePts[6]) || (!thesePts[2] && thesePts[6])){
						thisPoint = getBooleanPrimGeometry(p2a, p2b, p6a, p6b, hash1, primSet, c2, c6, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[6] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Z side, get the high Z edge
				if(k == zdim-1){
				///EDGE 2//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[1] && !thesePts[3]) || (!thesePts[1] && thesePts[3])){
						thisPoint = getBooleanPrimGeometry(p1a, p1b, p3a, p3b, hash1, primSet, c1, c3, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[2] = pointCounter; pointCounter++;
					}
				///EDGE 5//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[1] && !thesePts[5]) || (!thesePts[1] && thesePts[5])){
						thisPoint = getBooleanPrimGeometry(p1a, p1b, p5a, p5b, hash1, primSet, c1, c5, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[5] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high X and Y sides, get the high X and Y edge
				if(i == xdim-1 && j == ydim-1){
				///EDGE 11//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[6] && !thesePts[7]) || (!thesePts[6] && thesePts[7])){
						thisPoint = getBooleanPrimGeometry(p6a, p6b, p7a, p7b, hash1, primSet, c6, c7, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[11] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Y and Z sides, get the high Y and Z edge
				if(j == ydim-1 && k == zdim-1){
				///EDGE 7//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[3] && !thesePts[7]) || (!thesePts[3] && thesePts[7])){
						thisPoint = getBooleanPrimGeometry(p3a, p3b, p7a, p7b, hash1, primSet, c3, c7, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[7] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high X and Z sides, get the high X and Z edge
				if(k == zdim-1 && i == xdim-1){
				///EDGE 10//////////////////////////////////////////////////////////////////////////////////////
					if((thesePts[5] && !thesePts[7]) || (!thesePts[5] && thesePts[7])){
						thisPoint = getBooleanPrimGeometry(p5a, p5b, p7a, p7b, hash1, primSet, c5, c7, op);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[10] = pointCounter; pointCounter++;
					}
				}
								
//				if(oldPtCounter != pointCounter){ 
//					printf("NEW POINTS ADDED! oldPtCounter: %i pointCounter: %i\n", oldPtCounter, pointCounter); 
//				}

//				///RESET EARLIER EDGES BECAUSE WE DID NOT SET THEM BEFORE
				int * tempCube;
				if(theseEdges[0] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[8] = theseEdges[0];
					}
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[3] = theseEdges[0];
					}
					if(i != 0 && j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[11] = theseEdges[0];
					}
				}
				if(theseEdges[1] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[9] = theseEdges[1];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[2] = theseEdges[1];
					}
					if(i != 0 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[10] = theseEdges[1];
					}
				}
				if(theseEdges[2] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[10] = theseEdges[2];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[1] = theseEdges[2];
					}
					if(i != 0 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[9] = theseEdges[2];
					}
				}
				if(theseEdges[3] != -1){
					if(i != 0){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[11] = theseEdges[3];
					}
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[0] = theseEdges[3];
					}
					if(i != 0 && j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i-1)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[8] = theseEdges[3];
					}
				}
				if(theseEdges[4] != -1){
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[6] = theseEdges[4];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[5] = theseEdges[4];
					}
					if(j != 0 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k-1);
						tempCube[7] = theseEdges[4];
					}
				}
				if(theseEdges[5] != -1){
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[7] = theseEdges[5];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[4] = theseEdges[5];
					}
					if(j != 0 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+1);
						tempCube[6] = theseEdges[5];
					}
				}
				if(theseEdges[6] != -1){
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[4] = theseEdges[6];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[7] = theseEdges[6];
					}
					if(j != ydim-1 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k-1);
						tempCube[5] = theseEdges[6];
					}
				}
				if(theseEdges[7] != -1){
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[5] = theseEdges[7];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[6] = theseEdges[7];
					}
					if(j != ydim-1 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+1);
						tempCube[4] = theseEdges[7];
					}
				}
				if(theseEdges[8] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[0] = theseEdges[8];
					}
					if(j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[11] = theseEdges[8];
					}
					if(i != xdim-1 && j != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j-1)*zdim) + k+0);
						tempCube[3] = theseEdges[8];
					}
				}
				if(theseEdges[9] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[1] = theseEdges[9];
					}
					if(k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[10] = theseEdges[9];
					}
					if(i != xdim-1 && k != 0){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k-1);
						tempCube[2] = theseEdges[9];
					}
				}
				if(theseEdges[10] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[2] = theseEdges[10];
					}
					if(k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[9] = theseEdges[10];
					}
					if(i != xdim-1 && k != zdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+1);
						tempCube[1] = theseEdges[10];
					}
				}
				if(theseEdges[11] != -1){
					if(i != xdim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+0)*zdim) + k+0);
						tempCube[3] = theseEdges[11];
					}
					if(j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+0)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[8] = theseEdges[11];
					}
					if(i != xdim-1 && j != ydim-1){
						tempCube = (int*) mapsto_set(cubes, ((i+1)*ydim*zdim) + ((j+1)*zdim) + k+0);
						tempCube[0] = theseEdges[11];
					}
				}

				//delete[](dist);
				delete[](c0); delete[](c1); delete[](c2); delete[](c3);
				delete[](c4); delete[](c5); delete[](c6); delete[](c7);
				progress_cube++;
			}
		}		
	}

	///printf("POINTS GENERATED: %i\n", size_set(points));

	////Second Pass: generate Triangle List.
	printf("Generating triangles...\n"); fflush(stdout);
	int triListSize = 0;
	int triListCap = 100;
	int ** triList = new int*[triListCap];
	double ** normList = new double*[triListCap];
	vertTriSet * triSet = new vertTriSet();

	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				///look up the internal/external points, as stored before.
				int currentCubeIndex = (i*ydim*zdim) + (j*zdim) + k;
				bool * thesePts = (bool *) mapsto_set(corners, currentCubeIndex);
				int p0i = 0; if(thesePts[0]){ p0i = 1; }
				int p1i = 0; if(thesePts[1]){ p1i = 1; }
				int p2i = 0; if(thesePts[2]){ p2i = 1; }
				int p3i = 0; if(thesePts[3]){ p3i = 1; }
				int p4i = 0; if(thesePts[4]){ p4i = 1; }
				int p5i = 0; if(thesePts[5]){ p5i = 1; }
				int p6i = 0; if(thesePts[6]){ p6i = 1; }
				int p7i = 0; if(thesePts[7]){ p7i = 1; }

				///figure out what case we in the cubetable we will use.
				int cubeTableIndex = p0i + 2*p1i + 4*p2i + 8*p3i + 16*p4i + 32*p5i + 64*p6i + 128*p7i;
				int numPolysToAdd = table->sizes[cubeTableIndex];

				///create triangles for this cube.
				if(numPolysToAdd > 0){
					////expand triList and normList if they need to be expaded.
					if(triListSize + numPolysToAdd > triListCap){
						int ** newTriList = new int*[triListCap*2];
						double ** newNormList = new double*[triListCap*2];
						for(l = 0; l<triListSize; l++){
							newTriList[l] = triList[l];
							newNormList[l] = normList[l];
						}
						delete[](triList);
						delete[](normList);
						triList = newTriList;
						normList = newNormList;
						triListCap = triListCap*2;
					}
					///now create the triangles
//					printf("Creating %i Triangles.  CubeTableIndex: %i\n", numPolysToAdd, cubeTableIndex);
					int * polyList = table->polys[cubeTableIndex];
					int * thisCube = (int *) mapsto_set(cubes, currentCubeIndex);
					for(l = 0; l<numPolysToAdd; l++){
						int * newTriangle = new int[3];

						//printf("getting the point indices: numPolys: %i Edges: [%i %i %i]\n", numPolysToAdd, polyList[0], polyList[1], polyList[2] );
						///add the triangle into the triList.
						triList[triListSize] = newTriangle;

						///reversed order here to reverse norms - original polyList is reversed.
						newTriangle[0] = thisCube[polyList[3*l+0]];
						newTriangle[1] = thisCube[polyList[3*l+2]];	
						newTriangle[2] = thisCube[polyList[3*l+1]];

						//printf("adding indices to the triList\n");
						//store point indices into triList, for this triangle
						triSet->addEdge(newTriangle[0], triListSize);	
						triSet->addEdge(newTriangle[1], triListSize);
						triSet->addEdge(newTriangle[2], triListSize);

						//printf("Getting point coords: newTriangle[0]: %i, newTriangle[1]: %i newTriangle[2]: %i\n", newTriangle[0], newTriangle[1], newTriangle[2]);
						//get point Vals (do not delete these)
						double * firstPt = (double *) mapsto_set(points, newTriangle[0]);
						double * secndPt = (double *) mapsto_set(points, newTriangle[1]);
						double * thirdPt = (double *) mapsto_set(points, newTriangle[2]);

						///get normals (delete first, second)
						//printf("Calcing first diff for X-prod\n");
						double * first = new double[3];
						first[0] = (secndPt[0] - firstPt[0]);
						first[1] = (secndPt[1] - firstPt[1]);
						first[2] = (secndPt[2] - firstPt[2]);

						//printf("Calcing 2nd diff for X-prod\n");
						double * second = new double[3];
						second[0] = (thirdPt[0] - firstPt[0]);
						second[1] = (thirdPt[1] - firstPt[1]);
						second[2] = (thirdPt[2] - firstPt[2]);

						//printf("Calcing X-prod\n");
						double * normal = crossProd(first, second);
						double * normalizedNormal = normalizeVector(normal);
						normList[triListSize] = normalizedNormal;

						////clear and increment
						delete[](first);
						delete[](second);
						delete[](normal);

						triListSize++;
					}
				}
				else{
					///commented out bc it pritns out too much
					//printf("No Tris To Make!\n");
				}
				
			}
		}
	}

	////Third Pass: produce actual triangles and normals.
	printf("Generating Triangles and Normals\n"); fflush(stdout);
	int numPts;
	double * pts;
	double * norms;
	int numTris;
	int * tris;
	double * triangleNorms;

	double reverseNorm = 1;

	///Copy over the points	
	numPts = size_set(points);
	pts = new double[3*size_set(points)];
	for(i = 0; i<size_set(points); i++){
		double * tempPt = (double *) mapsto_set(points, i);
		pts[3*i+0] = tempPt[0];	pts[3*i+1] = tempPt[1];	pts[3*i+2] = tempPt[2];
	}
	///Copy over the triangles and triangle normals
	numTris = triListSize;
	tris = new int[3*triListSize];
	triangleNorms = new double[3*triListSize];
	for(i = 0; i<triListSize; i++){
		tris[3*i+0] = triList[i][0];	
		tris[3*i+1] = triList[i][1];	
		tris[3*i+2] = triList[i][2];
		triangleNorms[3*i+0] = reverseNorm*normList[i][0];	
		triangleNorms[3*i+1] = reverseNorm*normList[i][1];	
		triangleNorms[3*i+2] = reverseNorm*normList[i][2];
	}
	
	////STEP 3: Generate Averaged normals
	norms = new double[3*numPts];///this is output
	for(i = 0; i<triListSize; i++){
		int * thisTri = triList[i];
		////for each vertex in the triangle, get the list of adjacent triangles, and get the normals.
		///This loop processes each point in the triangle separately.
		for(j = 0; j<3; j++){
			int thisPt = thisTri[j];
//			double * thisVector = (double *) mapsto_set(points, thisPt);  ///this gets the coordinates of the point
			set_t neighbors = triSet->getNeighbors(thisPt); ///this gets the triangles that touch this point.
			int numNeigh = size_set(neighbors);  ///this is the number of triangles adjacent
			double * normal = new double[3];	///the normal for this point
			normal[0] = 0;		normal[1] = 0;		normal[2] = 0;
			for(k = 0; k<numNeigh; k++){
				normal[0] = normal[0] + normList[neighbors[k]][0];	///add up the normals of the adjacent traignels
				normal[1] = normal[1] + normList[neighbors[k]][1];
				normal[2] = normal[2] + normList[neighbors[k]][2];
			}
			normal[0] /= numNeigh;	///average the normals of adjacent triangles
			normal[1] /= numNeigh;
			normal[2] /= numNeigh;
			double * normedNormal = normalizeVector(normal);   //This is now the averaged normal; normalize it.

			norms[3*thisPt+0] = reverseNorm*normedNormal[0];	///Store the result into norms for this point in the triangle.
			norms[3*thisPt+1] = reverseNorm*normedNormal[1];
			norms[3*thisPt+2] = reverseNorm*normedNormal[2];

			delete[](normedNormal);
			delete[](normal);
		}
	}



	////STEP 4: Clean up all the memory allocations
	///////////////////////////////////////////////////////////////////////////////
	for(i = 0; i<size_set(points); i++){ double * vector = (double *) mapsto_set(points, i); delete[](vector); }
	for(i = 0; i<size_set(corners); i++){ int * corner = (int *) mapsto_set(corners, i); delete[](corner); }
	for(i = 0; i<size_set(cubes); i++){ int * cube = (int *) mapsto_set(cubes, i); delete[](cube); }
	for(i = 0; i<triListSize; i++){ delete[](triList[i]); delete[](normList[i]); }
	free_set(points);
	free_set(corners);
	free_set(cubes);
	delete[](triList);
	delete[](normList);
	delete(triSet);
	///////////////////////////////////////////////////////////////////////////////
	////STEP 4: Clean up all the memory allocations

	
	////STEP 5: Generate the output:	
	///////////////////////////////////////////////////////////////////////////////
	SurfaceObject * results;
	
	if(numPts == 0 || numTris == 0){ results = new SurfaceObject (); }
	results = new SurfaceObject (numPts, pts, norms, numTris, tris, triangleNorms);

	///debug
	//generateSURF(results, NULL, "TMESHTEST.SURF");
	///////////////////////////////////////////////////////////////////////////////
	////STEP 5: Generate the output:	
	

	////Step 6: Clear the output data:
	///////////////////////////////////////////////////////////////////////////////
	delete[](pts);
	delete[](norms);
	delete[](tris);
	delete[](triangleNorms);
	///////////////////////////////////////////////////////////////////////////////
	////Step 6: Clear the output data:


	/////Return the result
	///////////////////////////////////////////////////////////////////////////////
	return results;

}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



