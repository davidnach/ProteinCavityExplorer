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
 *       Local information for lattice hashing
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

#include "LatticeObj.h"


///constructor for latticeObj
///  cors fills in for corners
///  ind is the index of thsi latticeObj
///  sobj is the surfaceObj
LatticeObj::LatticeObj( double * cors, int ind, SurfaceObject * sobj )
{
	int i = 0;
	index = ind;

	///set the 8*3 corner coords
	corners = new double[24];
	for(i = 0; i<24; i++){ corners[i] = cors[i]; }

	indices = alloc_set(0);

	/// get the extents of the cube.
/*	maxX=-HUGE_VAL;	minX=HUGE_VAL;
	maxY=-HUGE_VAL;	minY=HUGE_VAL;
	maxZ=-HUGE_VAL;	minZ=HUGE_VAL;
	for(i = 0; i<8; i++){
		if(corners[3*i+0] > maxX){ maxX = corners[3*i+0]; }
		if(corners[3*i+0] < minX){ minX = corners[3*i+0]; }
		if(corners[3*i+1] > maxY){ maxY = corners[3*i+1]; }
		if(corners[3*i+1] < minY){ minY = corners[3*i+1]; }
		if(corners[3*i+2] > maxZ){ maxZ = corners[3*i+2]; }
		if(corners[3*i+2] < minZ){ minZ = corners[3*i+2]; }
	}
*/
	surf = sobj;
	adj = NULL;
	highlight = NOT_HIGHLIGHTED;
}


/*
LatticeObj::LatticeObj(  )
{
	corners = new double[0];
	indices = alloc_set(0);
	adj = alloc_set(0);
}
*/


LatticeObj::~LatticeObj()
{
	delete[](corners);
	free_set(adj);
	free_set(indices);
}


///set the adjacency data;
///  a is the array of adjancies (always 27), subbing -1 for non-adjacencies.
///  one is the reflexibe adjacency
void LatticeObj::setAdjacency(int * a)
{
	int i = 0;

	if(adj != NULL){ free_set(adj); }

	///set the 27 adjacency coords (one is the reflexive adjacency)
	adj = alloc_set(0);
	for(i = 0; i<27; i++){
		if(a[i] != -1){
			adj = put_set(adj, a[i]);
		}
	}
}


///set the triangle data;
void LatticeObj::addTriangle(int t)
{
	indices = put_set(indices, t);
}




///determine if a triangle intersects the volume of the cube. (tcoords is 9 doubles = 3*3)
/*
this code uses the following cube layout, with 3D axes to the right.
      7---6            ^
      |\  |\       \   |
      | 4---5       +y |
      | | | |        \ z+
      3-|-2 |         \|  
       \|  \|          0--x+-->    
        0---1
*/

bool LatticeObj::triangleIntersect( double * tcoords )
{
	int i = 0;

	///this gets tested in two cases.
	bool result = false;

	/// get the extents of the cube.
	double maxX=-HUGE_VAL;	double minX=HUGE_VAL;
	double maxY=-HUGE_VAL;	double minY=HUGE_VAL;
	double maxZ=-HUGE_VAL;	double minZ=HUGE_VAL;
	for(i = 0; i<8; i++){
		if(corners[3*i+0] > maxX){ maxX = corners[3*i+0]; }
		if(corners[3*i+0] < minX){ minX = corners[3*i+0]; }
		if(corners[3*i+1] > maxY){ maxY = corners[3*i+1]; }
		if(corners[3*i+1] < minY){ minY = corners[3*i+1]; }
		if(corners[3*i+2] > maxZ){ maxZ = corners[3*i+2]; }
		if(corners[3*i+2] < minZ){ minZ = corners[3*i+2]; }
	}

	result = triCubeIntersect(tcoords, minX, maxX, minY, maxY, minZ, maxZ);
	return result;
}


/*
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
///This old apparatus was developed to test a new implementation of segmentIntersect
///and is now deprecated.
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool LatticeObj::segmentIntersect( double * a, double * b )
{
//	bool test1 = segmentIntersectOld( a, b );
	bool test2 = segmentIntersectNew( a, b );
//	if(test1 != test2){
//		printf("WTF!!!! SEGMENT INTERSECTION DOES NOT MATCH!!!\n");
//		printf("a[0]=%f; a[1]=%f; a[2]=%f;  b[0]=%f; b[1]=%f; b[2]=%f;\n", a[0], a[1], a[2], b[0], b[1], b[2] );
//		printf("minX=%f; minY=%f; minZ=%f;   maxX=%f; maxY=%f; maxZ=%f; \n", 
//			corners[3*0+0], corners[3*0+1], corners[3*0+2], corners[3*6+0], corners[3*6+1], corners[3*6+2]);
//		printf("test1: %i  test2: %i\n", test1, test2);
//	}
	return test2;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
///This old apparatus was developed to test a new implementation of segmentIntersect
///and is now deprecated.
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
///determine if a segment intersects the volume of the cube.
bool LatticeObj::segmentIntersectOld( double * a, double * b )
{
	int i = 0;

	///this gets tested in two cases.
	bool result = false;

	/// get the extents of the cube.
	double maxX=-HUGE_VAL;	double minX=HUGE_VAL;
	double maxY=-HUGE_VAL;	double minY=HUGE_VAL;
	double maxZ=-HUGE_VAL;	double minZ=HUGE_VAL;
	for(i = 0; i<8; i++){
		if(corners[3*i+0] > maxX){ maxX = corners[3*i+0]; }
		if(corners[3*i+0] < minX){ minX = corners[3*i+0]; }
		if(corners[3*i+1] > maxY){ maxY = corners[3*i+1]; }
		if(corners[3*i+1] < minY){ minY = corners[3*i+1]; }
		if(corners[3*i+2] > maxZ){ maxZ = corners[3*i+2]; }
		if(corners[3*i+2] < minZ){ minZ = corners[3*i+2]; }
	}

	//First, test if any one of the segment endpoints is in the cube.
	///if any point is within the bounds, the triangle is within the cube.
	if( (a[0] > minX) && (a[0] < maxX) && (a[1] > minY) && (a[1] < maxY) && (a[2] > minZ) && (a[2] < maxZ) ){ return true; }
	if( (b[0] > minX) && (b[0] < maxX) && (b[1] > minY) && (b[1] < maxY) && (b[2] > minZ) && (b[2] < maxZ) ){ return true; }

	//Second, test if the segment intersects with any of the quads of the cube.
	//cc stands for "corner Coords".
	double ** cc = new double*[8];
	for(i = 0; i<8; i++){
		cc[i] = new double[3];
		cc[i][0] = corners[3*i+0];
		cc[i][1] = corners[3*i+1];
		cc[i][2] = corners[3*i+2];
	}

	if(!result){ result = segQuadOverlapTest( a, b, cc[0], cc[3], cc[2], cc[1] ); }
	if(!result){ result = segQuadOverlapTest( a, b, cc[0], cc[1], cc[5], cc[4] ); }
	if(!result){ result = segQuadOverlapTest( a, b, cc[1], cc[2], cc[6], cc[5] ); }
	if(!result){ result = segQuadOverlapTest( a, b, cc[2], cc[3], cc[7], cc[6] ); }
	if(!result){ result = segQuadOverlapTest( a, b, cc[3], cc[0], cc[4], cc[7] ); }
	if(!result){ result = segQuadOverlapTest( a, b, cc[4], cc[5], cc[6], cc[7] ); }

	//clean up
	for(i = 0; i<8; i++){ delete[]( cc[i]); }
	delete[](cc);

	return result;
}
*/


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///determine if a segment intersects the volume of the cube. 
///(faster version - less allocations, less use of segQuadOverlap, which is slow)
///NOTE - THIS ONLY WORKS WITH THE corners ARRAY LAID OUT AS IT IS IN TriangleLatticeHash.
bool LatticeObj::segmentIntersect( double * a, double * b )
{
	///////////////////////////////////////////////////////////////////////////////////////
	///NOTE - THIS ONLY WORKS WITH THE corners ARRAY LAID OUT AS IT IS IN TriangleLatticeHash.
	///get the extents of the cube.
	double maxX=corners[3*6+0];	double minX=corners[3*0+0];
	double maxY=corners[3*6+1];	double minY=corners[3*0+1];
	double maxZ=corners[3*6+2];	double minZ=corners[3*0+2];

	///////////////////////////////////////////////////////////////////////////////////////
	///important to set these values incase the dimension is invariant over T, but WITHIN the interval
	double XbotT=-HUGE_VAL, XtopT=HUGE_VAL, YbotT=-HUGE_VAL, YtopT=HUGE_VAL, ZbotT=-HUGE_VAL, ZtopT=HUGE_VAL;
	
	/////////////////////////////////////////////
	///test X overlap////////////////////////////
	if( b[0]==a[0]){								///if the X dimension is invariant over t
		if( (a[0]<minX) || (a[0]>maxX) ){ return false; }	///and the X dimension is not within the interval, return false;
	}
	else{
		XbotT = (minX-a[0])/(b[0]-a[0]);				///Solve for the lowest point along the line a-b that will be in X range
		XtopT = (maxX-a[0])/(b[0]-a[0]);				///Solve for the highest point along the line a-b that will be in X range
		if(XbotT>XtopT){ double tmp = XbotT; XbotT=XtopT; XtopT=tmp; }	//swap the values if the intervals are switched.
		if( (XbotT>1) || (XtopT<0) ){ return false; }	///the interval that t must be in is not in the segment.
	}
	if(XbotT<0){XbotT=0;}							///constrain the parametrization range to the segment.
	if(XtopT>1){XtopT=1;}

	//////////////////////////////////////////////
	///test Y overlap/////////////////////////////
	if( b[1]==a[1]){ 								///if the Y dimension is invariant over t
		if( (a[1]<minY) || (a[1]>maxY) ){ return false; } ///and the Y dimension is not within the interval, return false;
	}
	else{
		YbotT = (minY-a[1])/(b[1]-a[1]);				///Solve for the lowest point along the line a-b that will be in Y range
		YtopT = (maxY-a[1])/(b[1]-a[1]);				///Solve for the highest point along the line a-b that will be in Y range
		if(YbotT>YtopT){ double tmp = YbotT; YbotT=YtopT; YtopT=tmp; }	//swap the values if the intervals are switched.
		if( (YbotT>1) || (YtopT<0) ){ return false; }	///the interval that t must be in is not in the segment.
	}
	if(YbotT<0){YbotT=0;}							///constrain the parametrization range to the segment.
	if(YtopT>1){YtopT=1;}

	//////////////////////////////////////////////////////////////////////////////
	///this is moved earlier so we can quit if these two intervals do not overlap.
	if( (YtopT < XbotT) || (YbotT > XtopT) ){ 
		return false;
	}

	//////////////////////////////////////////////
	///test Z overlap/////////////////////////////
	if( b[2]==a[2]){  								///if the Z dimension is invariant over t
		if( (a[2]<minZ) || (a[2]>maxZ) ){ return false; }	///and the Z dimension is not within the interval, return false;
	}
	else{
		ZbotT = (minZ-a[2])/(b[2]-a[2]);				///Solve for the lowest point along the line a-b that will be in Z range
		ZtopT = (maxZ-a[2])/(b[2]-a[2]);				///Solve for the highest point along the line a-b that will be in Z range
		if(ZbotT>ZtopT){ double tmp = ZbotT; ZbotT=ZtopT; ZtopT=tmp; }	//swap the values if the intervals are switched.
		if( (ZbotT>1) || (ZtopT<0) ){ return false; }	///the interval that t must be in is not in the segment.
	}
	if(ZbotT<0){ZbotT=0;}							///constrain the parametrization range to the segment.
	if(ZtopT>1){ZtopT=1;}

	//////////////////////////////////////////////////////////////////////////////////////////////
	///make sure that the t-intervals over which the segment is in the range of the cube overlaps.
	///if they do not overlap, then it fails.  We have already tested for X-Y overlap.  now test X-Z and Y-Z
	///we will quit if it is false, so getting to a test means that the others are true.
	
	if( (ZtopT < XbotT) || (ZbotT > XtopT) ){ return false; }
	if( (ZtopT < YbotT) || (ZbotT > YtopT) ){ return false; }
	
	return true;
}



///Finds the smallest distance to all triangles in the cube to the point
///returns HUGE_VAL if cube is empty.
///p is a double[3]
double LatticeObj::getClosestDist( double * p, int &tIndex )
{
	int i = 0;
//	int j = 0;
	double result = HUGE_VAL;
	double tempVal = 0;
	double * t;	///this will house the triangle coords we use;

	///allocate arrays if we will really use them.
	if( size_set(indices) == 0){ return result; }

	///find the distance to all triangles;
	for(i = 0; i<size_set(indices); i++){
		t = surf->getTriangle( indices[i] );

		///this calls the alternate "distanceTo()" which uses the 3x3 triangle definition.
		tempVal = distanceTo(p, t);
		if(tempVal < result){ 
			result = tempVal; 
			tIndex = indices[i];
		}

		delete[](t);
	}

	return result;
}


///Finds the intersection Pts that intersects segment a-b
//output array is packed by first having the size, then 3xVectors
//Input arrays are 3x vectors
double * LatticeObj::getIntersectPts( double * a, double * b)
{
	int i = 0;
//	int j = 0;
//	double minDist = HUGE_VAL; 
//	double dist;
	set_t ptBuffer = alloc_set(SP_MAP);
	
//	printf("Triangles tested: [");			///debug instrumentation - DO NOT DELETE
	
	for(i = 0; i<size_set(indices); i++){
//		printf("%i", indices[i]);			///debug instrumentation - DO NOT DELETE
		
		double * t1 = new double[3];
			t1[0] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+0])+0];
			t1[1] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+0])+1];
			t1[2] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+0])+2];
		double * t2 = new double[3];
			t2[0] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+1])+0];
			t2[1] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+1])+1];
			t2[2] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+1])+2];
		double * t3 = new double[3];		
			t3[0] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+2])+0];
			t3[1] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+2])+1];
			t3[2] = surf->surfacePoints[3*(surf->triangles[3*indices[i]+2])+2];
		double * pt = intersect_Seg_Triangle( a, b, t1, t2, t3);

		if(pt!= NULL){
//			printf("(Y <%f, %f, %f>), ", pt[0], pt[1], pt[2]);				///debug instrumentation - DO NOT DELETE
			bool fail = false;
			double txmin = a[0];  if(txmin > b[0]){ txmin = b[0]; }
			double tymin = a[1];  if(tymin > b[1]){ tymin = b[1]; }
			double tzmin = a[2];  if(tzmin > b[2]){ tzmin = b[2]; }
			double txmax = a[0];  if(txmax < b[0]){ txmax = b[0]; }
			double tymax = a[1];  if(tymax < b[1]){ tymax = b[1]; }
			double tzmax = a[2];  if(tzmax < b[2]){ tzmax = b[2]; }

			if( pt[0]>txmax||pt[0]<txmin || pt[1]>tymax||pt[1]<tymin || pt[2]>tzmax||pt[2]<tzmin ){
		//		/////////////debug (check to see if the point is not crazy)
		//		printf("****ERROR**** [TRI # %i]: SEGMENT-TRIANGLE INTERSECTION OUT OF BOUNDS\n", indices[i]);
		//		printf("segA: (%f %f %f)    segB: (%f %f %f)   \"INTERSECTION\": %f %f %f\n", 
		//			a[0], a[1], a[2], b[0], b[1], b[2], pt[0], pt[1], pt[2] );
		//		printf("T1: (%f %f %f)    T2: (%f %f %f)    T3: (%f %f %f)\n",
		//			t1[0], t1[1], t1[2], t2[0], t2[1], t2[2], t3[0], t3[1], t3[2] );
		//		/////////////debug (check to see if the point is not crazy)

				fail = true;
			}
			/////////////debug
			
			if(!fail){ ptBuffer = associate_set(ptBuffer, size_set(ptBuffer), pt); }
			else{ delete[](pt); }
		}
//		else{							///debug instrumentation - DO NOT DELETE
//			printf("(N), ");				///debug instrumentation - DO NOT DELETE
//		}								///debug instrumentation - DO NOT DELETE
		delete[](t1);	delete[](t2);	delete[](t3);
	}
//	printf("] {Total: %i}\n", size_set(ptBuffer));	///debug instrumentation - DO NOT DELETE



	double * result = new double[ 1 + (3*size_set(ptBuffer)) ];
	result[0] = size_set(ptBuffer);
	for(i = 0; i<size_set(ptBuffer); i++){
		double * tmp = (double *) mapsto_set(ptBuffer, i);
		result[(3*i+0)+1] = tmp[0];
		result[(3*i+1)+1] = tmp[1];
		result[(3*i+2)+1] = tmp[2];
		delete[](tmp);
	}
	free_set(ptBuffer);

	return result;
	
}



///determine if a triangle intersects the volume of the cube. (tcoords is 9 doubles = 3*3)
/*
this code uses the following cube layout, with 3D axes to the right.
      7---6            ^
      |\  |\       \   |
      | 4---5       +y |
      | | | |        \ z+
      3-|-2 |         \|  
       \|  \|          0--x+-->    
        0---1
*/

///Generalized from class function: determines if a triangle intersects with an axis aligned cube.
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
	cc[0] = new double[3];	cc[0][0] = xn;	cc[0][1] = yn;	cc[0][2] = zn;	//point 0   this code uses the following cube layout
	cc[1] = new double[3];	cc[1][0] = xp;	cc[1][1] = yn;	cc[1][2] = zn;	//point 1         7---6            ^                                                           
	cc[2] = new double[3];	cc[2][0] = xp;	cc[2][1] = yp;	cc[2][2] = zn;	//point 2         |\  |\       \   |                                                           
	cc[3] = new double[3];	cc[3][0] = xn;	cc[3][1] = yp;	cc[3][2] = zn;	//point 3         | 4---5       +y |                                                           
	cc[4] = new double[3];	cc[4][0] = xn;	cc[4][1] = yn;	cc[4][2] = zp;	//point 4         | | | |        \ z+                                                          
	cc[5] = new double[3];	cc[5][0] = xp;	cc[5][1] = yn;	cc[5][2] = zp;	//point 5         3-|-2 |         \|                                                           
	cc[6] = new double[3];	cc[6][0] = xp;	cc[6][1] = yp;	cc[6][2] = zp;	//point 6          \|  \|          0--x+-->                                                    
	cc[7] = new double[3];	cc[7][0] = xn;	cc[7][1] = yp;	cc[7][2] = zp;	//point 7           0---1                                                                      

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





