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
 * File: TriangleLatticeHash.cpp
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

//
/////////////////////////////////////////////////////

#include "TriangleLatticeHash.h"


///############################################################################################
///###   ## ### ### ### ### ###    ###   ##   # ### #   # ## ### ####   ## ## #################
///## ### #  #  ## # ## ### ### ### # ### ## ##  ## ## ## ## ## # ## ### # ## #################
///###  ###     # ### # ### ### ### # ### ## ##   # ## ##    # ### ##  ###    #################
///#####  # # # #     # ### ###    ## ### ## ## #   ## ## ## #     ####  # ## #################
///## ### # ### # ### # ### ### ##### ### ## ## ##  ## ## ## # ### # ### # ## #################
///###   ## ### # ### #   #   # ######   ##   # ### ## ## ## # ### ##   ## ## #################
///############################################################################################

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///This lightweight class keeps a nonredundant list of 3d points.
///It tests for nonredundancy by linear search, so it cannot
///    handle very many points without becoming slow.
SmallPointHash::SmallPointHash()
{
	points = alloc_set(SP_MAP);

}


SmallPointHash::~SmallPointHash()
{
	int i = 0;
	for(i = 0; i<size_set(points); i++){
		double * thisPt = (double *) mapsto_set(points, i);
		delete[](thisPt);
	}

	free_set(points);
}



void SmallPointHash::addPoint(double x, double y, double z)
{
	bool info = testPoint(x, y, z);
	
	//if it's not in the hash, add it.
	if(!info){
		double * newPt = new double[3];
		newPt[0] = x;
		newPt[1] = y;
		newPt[2] = z;

		points = associate_set(points, size_set(points), newPt);
	}
}

///wrapper
void SmallPointHash::addPoint(double * pt)
{
	addPoint(pt[0], pt[1], pt[2]);
}


/// CAUTION - LINEAR TEST.
///return true if we have it in the hash, false otherwise.
bool SmallPointHash::testPoint(double x, double y, double z)
{
	bool result = false;
	
	int i = 0;
	for(i = 0; i<size_set(points); i++){
		double * thisPt = (double *) mapsto_set(points, i);
		if(x != thisPt[0]){ continue; }
		else if(y != thisPt[1]){ continue; }
		else if(z != thisPt[2]){ continue; }
		else{
			//in this case, we have it
			return true;
		}
	}

	return result;
}


int SmallPointHash::size()
{
	return size_set(points);
}


///This returns a reference.  Not a new object.
double * SmallPointHash::getPoint(int index)
{
	if( index < 0 || index > size_set(points) ){
		printf("ERROR: Requested point from SmallPointHash Object outside of bounds.  Exitting.\n");
		exit(1);
	}
	
	double * pt = (double *) mapsto_set(points, index);
	return pt;
}











///#########################################################################################
///#   #   ##   ##  ###  ## ###   # ####  ##   #   #   ##  ##   # ## ##  ###  ## ## ########
///## ## ## ## ## ## # ## # ### ### ### ## ## ### ### ## ## # ### ## # ## # ## # ## ########
///## ## ## ## ## ## # #### ###  ## ### ## ## ### ### ## ####  ##    # ## ## ###    ########
///## ##   ### ##    # #  # ### ### ###    ## ### ### ## #### ### ## #    ### ## ## ########
///## ## ## ## ## ## # ## # ### ### ### ## ## ### ### ## ## # ### ## # ## # ## # ## ########
///## ## ## #   # ## ##  ##   #   #   # ## ## ### ##   ##  ##   # ## # ## ##  ## ## ########
///#########################################################################################
///this code uses the following cube layout, with 3D axes to the right.
///      7---6            ^
///      |\  |\       \   |
///      | 4---5       +y |
///      | | | |        \ z+
///      3-|-2 |         \|  
///       \|  \|          0--x+-->    
///        0---1
TriangleLatticeHash::TriangleLatticeHash(SurfaceObject * sobj, double resolution, char * n)
{
	printf("Generating Triangle Lattice Hash [%s]\n", n);
	name = new char[strlen(n)+10];
	sprintf(name, n);

	int i = 0;	int j = 0;	int k = 0;

	res = resolution;
	xdim = 0;	ydim = 0;	zdim = 0;
	surf = sobj;
	isInsideOut = false;
	num_cache_hits = 0;
	num_local_cubes_tests = 0;
	num_long_ray_tests = 0;
	num_highlit_cube_tests = 0;
	
	isEmpty = false;
	if( (sobj->numPoints == 0) || (sobj->numTriangles == 0) ){
		isEmpty = true;
		return;
	}

	///create the cache
	inCacheSize = TRIHASHCACHESIZE;
	ioCacheVecs = new double[ 3*inCacheSize ];	for(i = 0; i<3*inCacheSize; i++){ ioCacheVecs[i] = HUGE_VAL; }
	ioCacheBools = new bool[ inCacheSize ];
	newCachePtr = 0;
	oldCachePtr = 0;

	//set initial extents
	xmax=-HUGE_VAL;	xmin=HUGE_VAL;
	ymax=-HUGE_VAL;	ymin=HUGE_VAL;
	zmax=-HUGE_VAL;	zmin=HUGE_VAL;
	//compute extents
	for(i = 0; i<sobj->numPoints; i++){
		if(surf->surfacePoints[3*i+0] > xmax){ xmax = surf->surfacePoints[3*i+0]; }
		if(surf->surfacePoints[3*i+0] < xmin){ xmin = surf->surfacePoints[3*i+0]; }
		if(surf->surfacePoints[3*i+1] > ymax){ ymax = surf->surfacePoints[3*i+1]; }
		if(surf->surfacePoints[3*i+1] < ymin){ ymin = surf->surfacePoints[3*i+1]; }
		if(surf->surfacePoints[3*i+2] > zmax){ zmax = surf->surfacePoints[3*i+2]; }
		if(surf->surfacePoints[3*i+2] < zmin){ zmin = surf->surfacePoints[3*i+2]; }
	}

	///automatically expand the size by LATTICE_PADDING in all directions;
	xmax = xmax + (LATTICE_PADDING*res);	xmin = xmin - (LATTICE_PADDING*res);
	ymax = ymax + (LATTICE_PADDING*res);	ymin = ymin - (LATTICE_PADDING*res);
	zmax = zmax + (LATTICE_PADDING*res);	zmin = zmin - (LATTICE_PADDING*res);

	///now compute the dimensions of the grid
	xdim = ( (int) ((xmax-xmin)/res)) + 1;
	ydim = ( (int) ((ymax-ymin)/res)) + 1;
	zdim = ( (int) ((zmax-zmin)/res)) + 1;

	///now update the max's to correlate to the dimension of the grid
	xmax = xmin + res*((int)xdim);
	ymax = ymin + res*((int)ydim);
	zmax = zmin + res*((int)zdim);

	///now generate the grid.
	objs = new LatticeObj*[xdim*ydim*zdim];
	int counter = 0;
	grid = new int**[xdim];
	for(i = 0; i<xdim; i++){
		grid[i] = new int*[ydim];
		for(j = 0; j<ydim; j++){
			grid[i][j] = new int[zdim];
			for(k = 0; k<zdim; k++){
				double * cors = new double[8*3];	//corner coordinates
				double xbot = xmin+(i+0)*res;		double xtop = xmin+(i+1)*res;
				double ybot = ymin+(j+0)*res;		double ytop = ymin+(j+1)*res;
				double zbot = zmin+(k+0)*res;		double ztop = zmin+(k+1)*res;

			//	cube 0				cube 1				cube 2				cube 3
				cors[3*0+0] = xbot;	cors[3*1+0] = xtop;	cors[3*2+0] = xtop;	cors[3*3+0] = xbot;
				cors[3*0+1] = ybot;	cors[3*1+1] = ybot;	cors[3*2+1] = ytop;	cors[3*3+1] = ytop;
				cors[3*0+2] = zbot;	cors[3*1+2] = zbot;	cors[3*2+2] = zbot;	cors[3*3+2] = zbot;

			//	cube 4				cube 5				cube 6				cube 7
				cors[3*4+0] = xbot;	cors[3*5+0] = xtop;	cors[3*6+0] = xtop;	cors[3*7+0] = xbot;
				cors[3*4+1] = ybot;	cors[3*5+1] = ybot;	cors[3*6+1] = ytop;	cors[3*7+1] = ytop;
				cors[3*4+2] = ztop;	cors[3*5+2] = ztop;	cors[3*6+2] = ztop;	cors[3*7+2] = ztop;

				objs[counter] = new LatticeObj( cors, counter, surf );
				grid[i][j][k] = counter;
				delete[](cors);
				counter++;
			}
		}
	}

	build();
}




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
TriangleLatticeHash::~TriangleLatticeHash()
{
	int i = 0;
	int j = 0;
	
	delete[](name);
	delete[](ioCacheVecs);
	delete[](ioCacheBools);
	
	for(i = 0; i<xdim*ydim*zdim; i++){
		delete(objs[i]);
	}
	delete[](objs);

	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			delete[](grid[i][j]);
		}
		delete[](grid[i]);
	}
	delete[](grid);

}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///sets the cache size to a specified size.
void TriangleLatticeHash::setCacheSize(int s)
{
	inCacheSize = s;
	delete[](ioCacheVecs);
	ioCacheVecs = new double[ 3*inCacheSize ];
	delete[](ioCacheBools);
	ioCacheBools = new bool[ inCacheSize ];
	newCachePtr = 0;
	oldCachePtr = 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///if this is out of bounds, return -1;
///if this is not out of bounds, return the right int.
int TriangleLatticeHash::getCube(int x, int y, int z)
{
	if( (x<0) || (x>xdim-1) ){ return -1; }	
	if( (y<0) || (y>ydim-1) ){ return -1; }	
	if( (z<0) || (z>zdim-1) ){ return -1; }	

	return grid[x][y][z];
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///returns the cubes which this point intersects
bool TriangleLatticeHash::inRange(int dim, int val)
{
	bool result = true;
	if(dim == 0){
		if( (val < 0) || (val > xdim-1) ){ result = false; }
	}
	if(dim == 1){
		if( (val < 0) || (val > ydim-1) ){ result = false; }
	}
	if(dim == 2){
		if( (val < 0) || (val > zdim-1) ){ result = false; }
	}
	return result;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///returns the cubes which this triangle intersects.
set_t TriangleLatticeHash::getCubes(double * t)
{
	set_t result = alloc_set(0);
	int i = 0;	int j = 0;	int k = 0;
	
	/////Find the triangle bounds
	double sxmax = -HUGE_VAL;	double sxmin = HUGE_VAL;
	double symax = -HUGE_VAL;	double symin = HUGE_VAL;
	double szmax = -HUGE_VAL;	double szmin = HUGE_VAL;
	for(i = 0; i<3; i++){
		if(sxmax < t[3*i+0]){ sxmax = t[3*i+0]; }
		if(sxmin > t[3*i+0]){ sxmin = t[3*i+0]; }
		if(symax < t[3*i+1]){ symax = t[3*i+1]; }
		if(symin > t[3*i+1]){ symin = t[3*i+1]; }
		if(szmax < t[3*i+2]){ szmax = t[3*i+2]; }
		if(szmin > t[3*i+2]){ szmin = t[3*i+2]; }
	}

	/////map the indices as if the grid was infinite.
	int ihi, ilo, jhi, jlo, khi, klo;
	ihi = (int)( (sxmax-xmin)/res );
	jhi = (int)( (symax-ymin)/res );
	khi = (int)( (szmax-zmin)/res );
	if( fmod( (sxmin-xmin), res ) == 0 ){ ilo = ((int)((sxmin-xmin)/res)) -1; }
	else{ ilo = ((int)((sxmin-xmin)/res)); }
	if( fmod( (symin-ymin), res ) == 0 ){ jlo = ((int)((symin-ymin)/res)) -1; }
	else{ jlo = ((int)((symin-ymin)/res)); }
	if( fmod( (szmin-zmin), res ) == 0 ){ klo = ((int)((szmin-zmin)/res)) -1; }
	else{ klo = ((int)((szmin-zmin)/res)); }

	///bound the indices;  These are inclusive indices;
	if( ihi >= xdim ){ ihi = xdim-1; }
	if( ilo < 0 ){ ilo = 0; }
	if( jhi >= ydim ){ jhi = ydim-1; }
	if( jlo < 0 ){ jlo = 0; }
	if( khi >= zdim ){ khi = zdim-1; }
	if( klo < 0 ){ klo = 0; }

	///retrieve the cubes
	for(i = ilo; i<=ihi; i++){
		for(j = jlo; j<=jhi; j++){
			for(k = klo; k<=khi; k++){
				result = put_set(result, grid[i][j][k] );
			}
		}
	}

	return result;
}





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///build the hash
///now we add the adjacency info and the triangles
void TriangleLatticeHash::build()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;
	int n = 0;

	int * a = new int[27];	//adjacency data
//	int numCubes = xdim*ydim*zdim;
	
	///iterate through cubes, generate adjacency data. (iterate through the cubes in 3 dimensions)
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				///process neighbors (in 3 dimensions)
				for(l = -1; l<2; l++){
					for(m = -1; m<2; m++){
						for(n = -1; n<2; n++){
							///getCube inserts -1 when there is no neighbor in some location
							a[ (l+1)*9 + (m+1)*3 + (n+1) ] = getCube(i+l, j+m, k+n);
						}
					}
				}
				///set the adjacency and triangles.
				LatticeObj * tempObj = objs[i*(ydim*zdim) + j*zdim + k];
				///setAdjacency knows that -1 means there is nothing there, and doesnt try to add -1.
				tempObj->setAdjacency(a);
			}
		}
	}
	delete[](a);

	///iterate through triangles, map into cubes.
	double * t = new double[9];
	for(i = 0; i<surf->numTriangles; i++){
		t[0] = surf->surfacePoints[ 3*surf->triangles[3*i+0]+0 ];
		t[1] = surf->surfacePoints[ 3*surf->triangles[3*i+0]+1 ];
		t[2] = surf->surfacePoints[ 3*surf->triangles[3*i+0]+2 ];
		t[3] = surf->surfacePoints[ 3*surf->triangles[3*i+1]+0 ];
		t[4] = surf->surfacePoints[ 3*surf->triangles[3*i+1]+1 ];
		t[5] = surf->surfacePoints[ 3*surf->triangles[3*i+1]+2 ];
		t[6] = surf->surfacePoints[ 3*surf->triangles[3*i+2]+0 ];
		t[7] = surf->surfacePoints[ 3*surf->triangles[3*i+2]+1 ];
		t[8] = surf->surfacePoints[ 3*surf->triangles[3*i+2]+2 ];

		set_t myCubes = getCubes(t);
		for(j = 0; j<size_set(myCubes); j++){
			LatticeObj * tempObj = objs[myCubes[j]];
			if( tempObj->triangleIntersect( t ) ){
				tempObj->addTriangle(i);
			}
		}
		free_set(myCubes);
	}

	delete[](t);
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///this finds the distance to the closest triangle (this is a query)
/// returns HUGE_VAL if there are no triangles.
///RECODED with sphereical corner checking.  adjacency stuff is now deprecated.
set_t TriangleLatticeHash::distTo(double * pt, double * outputDist)
{
	int i = 0;
	int j = 0;

	set_t result = alloc_set(0);
	if(isEmpty){ return result; }
	
	double radius = res;
	
	///rename the points
	double x = pt[0];	double y = pt[1];	double z = pt[2];

	////out of bounds error checking - guarantees that the radius hits at least one corner of one cube.
	bool outOfBounds = false;	
	bool xlo = x<xmin;	bool xhi = x>xmax;
	bool ylo = y<ymin;	bool yhi = y>ymax;
	bool zlo = z<zmin;	bool zhi = z>zmax;
	if( xlo || xhi || ylo || yhi || zlo || zhi ){ outOfBounds = true; }
	if(outOfBounds){
		double xdist = 0.0;		if(xlo){ xdist = xmin-x; }		if(xhi){ xdist = x-xmax; }
		double ydist = 0.0;		if(ylo){ ydist = ymin-y; }		if(yhi){ ydist = y-ymax; }
		double zdist = 0.0;		if(zlo){ zdist = zmin-z; }		if(zhi){ zdist = z-zmax; }
		double boundsDist = sqrt( xdist*xdist + ydist*ydist + zdist*zdist );
		if(boundsDist > radius){ radius = boundsDist+res; }
	}

	///expand the testing if we need to.
	///Do not count degenerate triangles
	set_t testCubes = NULL;
	set_t testTriangles = alloc_set(0);
	while( (testCubes == NULL) || (size_set(testTriangles) == 0)){
		///get the cubes in the current radius
		testCubes = getCubes(pt, radius);
		///for each cube, count the number of nondegenerate triangles
		for(i = 0; i<size_set( testCubes ); i++){
			for(j = 0; j<size_set(objs[testCubes[i]]->indices); j++){
				///These tests are now taken care of by the functions that call this method
				///This ensures that methods that require accuracy for distance get accuracy,
				/// and that methods that require the nearest triangles for normal calculations (isInside())
				/// can filter the output as they please.
				int triangle_A = surf->triangles[3*objs[testCubes[i]]->indices[j]+0];
				int triangle_B = surf->triangles[3*objs[testCubes[i]]->indices[j]+1];
				int triangle_C = surf->triangles[3*objs[testCubes[i]]->indices[j]+2];
				double thisArea = triangleArea(surf->surfacePoints, triangle_A, triangle_B, triangle_C);
			//	bool colinear = isColinear(surf->surfacePoints, triangle_A, triangle_B, triangle_C);
			//	if( thisArea > 0 && !colinear ){ 
				if( thisArea > 0 ){ 
					testTriangles = put_set(testTriangles, objs[testCubes[i]]->indices[j]);
				}
			}
		}
		if( size_set(testTriangles) == 0 ){ 
			free_set(testCubes);
			//this added in as a precaution. BYC 05-07-2012
			testCubes = NULL;
		}
		radius += res;
	}
	free_set(testCubes);

	if(size_set(testTriangles)==0){
		printf("WTF!!!\n");
	}

	///find teh closest Triangle
	///Do not count degenerate triangles - No degenerate triangles are allowed to this stage
	///using the step above.
	double * ranges = new double[size_set(testTriangles)];
	int * refs = new int[size_set(testTriangles)];
	for(i = 0; i<size_set(testTriangles); i++){
		double * tri = surf->getTriangle( testTriangles[i] );
		ranges[i] = distanceTo( pt, tri );
		refs[i] = testTriangles[i];
		delete[](tri);
	}
	mergeSort(ranges, refs, size_set(testTriangles));
	double closest = ranges[0];
	i = 0;
	while( (i<size_set(testTriangles)) && (ranges[i] < closest+.00001)){
		result = put_set(result, refs[i]);
		i++;
	}
	delete[](ranges);
	delete[](refs);
	free_set(testTriangles);

	//printf("Closest Dist: %f\n", closest);
	if(outputDist != NULL){
		outputDist[0] = closest;
	}
	return result;
}




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///this code uses the following cube layout, with 3D axes to the right.
///      7---6            ^
///      |\  |\       \   |
///      | 4---5       +y |
///      | | | |        \ z+
///      3-|-2 |         \|  
///       \|  \|          0--x+-->    
///        0---1
///
///returns the cubes which this sphere intersects
set_t TriangleLatticeHash::getCubes(double * pt, double rad)
{
	int i, j, k;
	set_t result = alloc_set(0);
	if(isEmpty){ return result; }

	///rename the points
	double x = pt[0];	double y = pt[1];	double z = pt[2];

	////out of bounds error checking
//	bool outOfBounds = false;	
	bool xlo = x<xmin;
	bool xhi = x>xmax;
	bool ylo = y<ymin;
	bool yhi = y>ymax;
	bool zlo = z<zmin;
	bool zhi = z>zmax;
	///if we are out of bounds, check the radius.
	if( xlo || xhi || ylo || yhi || zlo || zhi ){ 
		double * p0 = new double[3];double * p1 = new double[3];double * p2 = new double[3];double * p3 = new double[3];
		p0[0] = xmin;	p1[0] = xmax;	p2[0] = xmax;	p3[0] = xmin;
		p0[1] = ymin;	p1[1] = ymin;	p2[1] = ymax;	p3[1] = ymax;
		p0[2] = zmin;	p1[2] = zmin;	p2[2] = zmin;	p3[2] = zmin;

		double * p4 = new double[3];double * p5 = new double[3];double * p6 = new double[3];double * p7 = new double[3];
		p4[0] = xmin;	p5[0] = xmax;	p6[0] = xmax;	p7[0] = xmin;
		p4[1] = ymin;	p5[1] = ymin;	p6[1] = ymax;	p7[1] = ymax;
		p4[2] = zmax;	p5[2] = zmax;	p6[2] = zmax;	p7[2] = zmax;

		double range = HUGE_VAL; double tempRange;
		tempRange = distanceTo(pt, p0, p1, p2 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p0, p2, p3 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p0, p4, p7 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p0, p7, p3 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p0, p1, p5 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p0, p5, p4 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p1, p2, p6 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p1, p6, p5 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p2, p3, p7 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p2, p7, p6 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p4, p7, p6 );	if(tempRange < range){ range = tempRange; }
		tempRange = distanceTo(pt, p4, p6, p5 );	if(tempRange < range){ range = tempRange; }

		delete[](p0);delete[](p1);delete[](p2);delete[](p3);
		delete[](p4);delete[](p5);delete[](p6);delete[](p7);

		if(range > rad){
			return result;
		}
	}

	/////Find the sphere bounds
	double sxmax = x+rad;
	double sxmin = x-rad;
	double symax = y+rad;
	double symin = y-rad;
	double szmax = z+rad;
	double szmin = z-rad;

	/////map the indices as if the grid was infinite.
	int ihi, ilo, jhi, jlo, khi, klo;
	ihi = (int)( (sxmax-xmin)/res );
	jhi = (int)( (symax-ymin)/res );
	khi = (int)( (szmax-zmin)/res );
	if( fmod( (sxmin-xmin), res ) == 0 ){ ilo = ((int)((sxmin-xmin)/res)) -1; }
	else{ ilo = ((int)((sxmin-xmin)/res)); }
	if( fmod( (symin-ymin), res ) == 0 ){ jlo = ((int)((symin-ymin)/res)) -1; }
	else{ jlo = ((int)((symin-ymin)/res)); }
	if( fmod( (szmin-zmin), res ) == 0 ){ klo = ((int)((szmin-zmin)/res)) -1; }
	else{ klo = ((int)((szmin-zmin)/res)); }

	///bound the indices;  These are inclusive indices;
	if( ihi >= xdim ){ ihi = xdim-1; }
	if( ilo < 0 ){ ilo = 0; }
	if( jhi >= ydim ){ jhi = ydim-1; }
	if( jlo < 0 ){ jlo = 0; }
	if( khi >= zdim ){ khi = zdim-1; }
	if( klo < 0 ){ klo = 0; }

	set_t tempCubes = alloc_set(0);
	for(i = ilo; i<=ihi; i++){
		for(j = jlo; j<=jhi; j++){
			for(k = klo; k<=khi; k++){
				tempCubes = put_set(tempCubes, grid[i][j][k] );
			}
		}
	}

	///now test cubes for penetration with the sphere.
	double * cornerPt = new double[3];
	result = alloc_set(0);
//	int debugCounter = 0;
	
//	printf("ilo: %i ihi: %i jlo: %i jhi: %i klo: %i khi: %i\n", ilo, ihi, jlo, jhi, klo, khi);

	///Note this is not a very complete test - faces and edges are not tested correctly here.
	for(i = 0; i<size_set(tempCubes); i++){
		LatticeObj * tempObj = objs[tempCubes[i]];
		for(j = 0; j<8; j++){
			cornerPt[0] = tempObj->corners[3*j+0];
			cornerPt[1] = tempObj->corners[3*j+1];
			cornerPt[2] = tempObj->corners[3*j+2];
			if( vectorSize( cornerPt[0]-x, cornerPt[1]-y, cornerPt[2]-z) < rad ){
				result = put_set(result, tempCubes[i]);
				break;
			}
		}
	}

	delete[](cornerPt);
	free_set(tempCubes);

	return result;
}


/*

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/// Inside Outside testing that uses Winding Number computation
double TriangleLatticeHash::isInsideWinding(double * pt)
{
	int i = 0;
	double result = 0;
	
	for(i = 0; i<surf->numTriangles; i++){
		double * t = surf->getTriangle(i); 
		double sphericalArea = sphereicalTriangleArea(pt, t);
		delete[](t);

		result += sphericalArea;
	}

//	////debug
//	if(result > 12){
//		result = 0;
//		for(i = 0; i<surf->numTriangles; i++){
//			double * t = surf->getTriangle(i); 
//			double sphericalArea = sphereicalTriangleArea(pt, t);
//			delete[](t);
//			printf("result += %f\n", sphericalArea);
//			result += sphericalArea;
//		}
//		
//		printf("TOTAL :%f\n", result);
//		exit(1);
//	}
//	////debug

	return result;
}

*/



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///this determines if the given point is inside the surface
///this is accomplished by casting a long ray (Segment) from pt, and counting the number
/// of logical intersections (identical intersection points are counted once).  Even
/// numbers of intersections mean that pt is outside, odd means that pt is inside.
///
//bool TriangleLatticeHash::isInsideRayTest(double * pt, double * dist)
//note that dist can be NULL
//
bool TriangleLatticeHash::isInsideRayTest(double * pt, double * dist)
{
	double * endPt = new double[3];

	///randomized ray:
	double endPta = (((double)(random()%2000)) - 1000.0)/1000.0;
	double endPtb = (((double)(random()%2000)) - 1000.0)/1000.0;
	double endPtc = (((double)(random()%2000)) - 1000.0)/1000.0;

	double * tempDir = normalizeVector(endPta, endPtb, endPtc);
	endPt[0] = pt[0] + tempDir[0]*1000;
	endPt[1] = pt[1] + tempDir[1]*1000;
	endPt[2] = pt[2] + tempDir[2]*1000;
	delete[](tempDir);
	
	//printf("pt[0]: %f    pt[1]: %f    pt[2]: %f\n", pt[0], pt[1], pt[2] );

	/*	
	///identify the closest direction to the border of hashed space
	double closestWall = HUGE_VAL;
	//0 xpos, 1 xneg, 2ypos, 3yneg, 4zpos, 5zneg
	int wallId = -1;
	if(xmax - pt[0] < closestWall){ wallId = 0; closestWall = xmax - pt[0];}
	if(pt[0] - xmin < closestWall){ wallId = 1; closestWall = pt[0] - xmin;}
	if(ymax - pt[1] < closestWall){ wallId = 2; closestWall = ymax - pt[1];}
	if(pt[1] - ymin < closestWall){ wallId = 3; closestWall = pt[1] - ymin;}
	if(zmax - pt[2] < closestWall){ wallId = 4; closestWall = zmax - pt[2];}
	if(pt[2] - zmin < closestWall){ wallId = 5; closestWall = pt[2] - zmin;}
	
	endPt[0] = pt[0]; endPt[1] = pt[1]; endPt[2] = pt[2]; 
	switch(wallId){
		case 0: endPt[0] = xmax; break;
		case 1: endPt[0] = xmin; break;
		case 2: endPt[1] = ymax; break;
		case 3: endPt[1] = ymin; break;
		case 4: endPt[2] = zmax; break;
		case 5: endPt[2] = zmin; break;
	}
	*/

	//	printf("a[0]=%f; a[1]=%f; a[2]=%f;  b[0]=%f; b[1]=%f; b[2]=%f;\n", pt[0], pt[1], pt[2], endPt[0], endPt[1], endPt[2] );

	int numInts = countIntersectionPts(pt, endPt, dist);

	//	printf("RAY TEST: num intersection pts: %i\n", numInts);

	double result = true;
	if( numInts%2 == 0){ result = false; }

	delete[](endPt);

	///get the most accurate distance;
	if(dist != NULL){
		set_t closeTris = distTo( pt, dist);
		free_set(closeTris);
	}

	return result;
}




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
////this wraps isInsideRayTest by testing several times, and deciding on the majority
////result, since we cannot depend on ray testing,
bool TriangleLatticeHash::isInsideRayTestMajority(double * pt, int numTests, bool * unanimous)
{
	int i = 0;
	int numInside = 0;
	for(i = 0; i<numTests; i++){
		bool val = isInsideRayTest(pt, NULL);
		if(val){ numInside++; }
	}

	bool result = false;
	if( ((double)numInside) > (((double)numTests)/2.0) ){
		//printf("Ray Test Inside: %i vs %i\n", numInside, numTests);
		result = true;
	}
//	else{ printf("Ray Test Outside: %i vs %i\n", numTests-numInside, numTests); }
	
	if( (numInside == 0) || (numInside==numTests) ){ unanimous[0] = true; }
	else{ unanimous[0] = false; }
	
	return result;
}




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///examines the closest triangle.  if the normal of that triangle points towards the point
/// then the point is outside (-1).  If it points away from the point, it is inside (1).  Otherwise
/// the test is inconclusive (0).  Results in parens.  If more than one triangle is closest,
/// the test is inconclusive.
int TriangleLatticeHash::isInsideClosestTriangle(double * pt)
{
	int result = 0;
	
	set_t closeTris = distTo( pt, NULL );
	if(size_set(closeTris) == 1){
		double * triangle = surf->getTriangle(closeTris[0]);
		
		double vector0 = triangle[0] - pt[0];
		double vector1 = triangle[1] - pt[1];
		double vector2 = triangle[2] - pt[2];

		double cross0 = ((triangle[4]-triangle[1])*(triangle[8]-triangle[2])) - ((triangle[5]-triangle[2])*(triangle[7]-triangle[1]));
		double cross1 = -( ((triangle[3]-triangle[0])*(triangle[8]-triangle[2])) - ((triangle[5]-triangle[2])*(triangle[6]-triangle[0])) );
		double cross2 = ((triangle[3]-triangle[0])*(triangle[7]-triangle[1])) - ((triangle[4]-triangle[1])*(triangle[6]-triangle[0]));

		double test = (vector0*cross0) + (vector1*cross1) + (vector2*cross2);
		if(fabs(test) < SMALL_NUM){ test = 0.0; }

		if(test>0){ result = 1; }
		if(test<0){ result = -1; }
	}	

	free_set(closeTris);
	return result;
}




///################################################################################################
///#     ##   ##     # ### ##   ##     #    ##     ## ## ##########################################
///### ### ### ### ###  ## # ### ### ### ### # ##### #### #########################################
///### ####  ##### ###   # ##  ##### ### ### #    ## #### #########################################
///### ######  ### ### #   ####  ### ### ### # ##### #### #########################################
///### ### ### ### ### ##  # ### ### ### ### # ##### #### #########################################
///#     ##   ##     # ### ##   ##     #    ##     ## ## ##########################################
///################################################################################################
///this is a convenient wrapper for the wrapper.
bool TriangleLatticeHash::isInside(double ptx, double pty, double ptz)
{
	double * input = new double[3];
	input[0] = ptx;
	input[1] = pty;
	input[2] = ptz;
	
	bool result = isInside(input);
	delete[](input);
	return result;
}


///this is the wrapper for insideOutside testing.
bool TriangleLatticeHash::isInside(double * pt)
{

	if(isEmpty){ 
		if(isInsideOut){ return true; }
		else{return false; }
	}

	//	printf("testing: %+lf %+lf %+lf\n", pt[0], pt[1], pt[2] );
	bool result;

	//check for a cached result;
	int val = isInsideCached(pt);
	if(val != 0){
		num_cache_hits++;	///benchmarking
		if(val == 1){return true;}
		else{return false;}
	}
	else{
		
		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		///Test here for out-of-bounds results
//		int boundsCheck = isinsideBoundsCheck(pt);
//		if( boundsCheck == -1 ){ return false; }
//		if( boundsCheck ==  1 ){ return true; }
		///Test here for out-of-bounds results
		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		///Test Here with segments to local cubes.
		int localCubesTest = isInsideLocalCubes(pt);
		///This func returns -1 if outside, +1 if inside, 0 if failure.
		if(localCubesTest != 0){
			num_local_cubes_tests++;	///benchmarking
			if(localCubesTest == -1){ result = false; }
			else{ result = true; }
		}
		///Test Here with segments to local cubes.
		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		////if local cubes fails, perform exaustive ray testing.
		else{
			num_long_ray_tests++;	///benchmarking
			bool unanimous;
			result = isInsideRayTestMajority(pt, NUM_INITIAL_RAY_TESTS, &unanimous);
			if(!unanimous){
				result = isInsideRayTestMajority(pt, NUM_THOROUGH_RAY_TESTS, &unanimous);
			}
		}
		////if local cubes fails, perform exaustive ray testing.
		///////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////

		///deposit in the cache
		ptCacheDeposit(pt, result);
	}	

	///flip the result if it is insideOut					///added to implement negative spaces
	if(isInsideOut){ result = !result; }

	return result;
}
///################################################################################################
///#     ##   ##     # ### ##   ##     #    ##     ## ## ##########################################
///### ### ### ### ###  ## # ### ### ### ### # ##### #### #########################################
///### ####  ##### ###   # ##  ##### ### ### #    ## #### #########################################
///### ######  ### ### #   ####  ### ### ### # ##### #### #########################################
///### ### ### ### ### ##  # ### ### ### ### # ##### #### #########################################
///#     ##   ##     # ### ##   ##     #    ##     ## ## ##########################################
///################################################################################################

///This returns 0 if inside the bounds, otherwise, -1 if not inside out (e.g. outside), 1 if inside out (e.g. inside).
/*
int TriangleLatticeHash::isinsideBoundsCheck(double * pt)
{
	int result = 0;
	
	if( pt[0] < xmin || pt[0] > xmax || pt[1] < ymin || pt[1] > ymax || pt[2] < zmin || pt[2] > zmax){
		if( isInsideOut ){
			result = 1;
		}
		else{
			result = -1;
		}
	}

	return result;
}
*/


int TriangleLatticeHash::isInsideCached(double * pt)
{
	int i = 0;
	bool oldResult;
	bool foundit = false;
	for(i = oldCachePtr; (i%inCacheSize)!=newCachePtr; i++){
		double px = ioCacheVecs[3*(i%inCacheSize)+0];
		double py = ioCacheVecs[3*(i%inCacheSize)+1];
		double pz = ioCacheVecs[3*(i%inCacheSize)+2];
		if( (pt[0]==px) && (pt[1]==py) && (pt[2]==pz) ){
			oldResult = ioCacheBools[(i%inCacheSize)];
			foundit = true;
			break;
		}
	}
	if(foundit){
	//	printf("cache hit!\n");
		if(oldResult){ return 1; }
		else{ return -1; }
	}
	else{
	//	printf("cache miss!\n");
		return 0;
	}
}

void TriangleLatticeHash::ptCacheDeposit(double * pt, bool ioVal)
{
	ioCacheVecs[3*newCachePtr+0] = pt[0];
	ioCacheVecs[3*newCachePtr+1] = pt[1];
	ioCacheVecs[3*newCachePtr+2] = pt[2];
	ioCacheBools[newCachePtr] = ioVal;

	newCachePtr++;
	if(newCachePtr == inCacheSize){ newCachePtr = 0; }
	if(newCachePtr == oldCachePtr){ oldCachePtr++;  }
	if(oldCachePtr == inCacheSize){ oldCachePtr = 0; }
}






//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
///################################################################################################///
///# ####   ###   #### ### ####   ## ### #    ##    ##   ### ## ###################################///
///# ### ### # ### ## # ## ### ### # ### # ### # #### ### # #### ##################################///
///# ### ### # ##### ### # ### ##### ### #    ##   ###  ### #### ##################################///
///# ### ### # #####     # ### ##### ### # ### # #######  # #### ##################################///
///# ### ### # ### # ### # ### ### # ### # ### # #### ### # #### ##################################///
///#   ##   ###   ## ### #   ##   ###   ##    ##    ##   ### ## ###################################///
///################################################################################################///
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
///This function tests inside/outside based on segment testing to nearby cubes.
///The only cubes tested are those that are highlighted as inside or outside.
///This func returns -1 if outside, +1 if inside, 0 if failure.
///  Failure happens two ways.  First, if there are no nearby empty cubes, we fail.
///		Second, if there is a contradiction in segment testing, we fial.
///		
///This algorithm assumes that the grid is mostly empty, that a hilighted cube is close by.		
///This will fail a lot if this assumption is not very true.
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
int TriangleLatticeHash::isInsideLocalCubes(double * pt)
{
	int i = 0;
	int j = 0;
	int result = 0;
	
	///we will check inside/outside status based on segment tests to nearby cubes.
	///LOCAL_CUBE_CHECK_MULTIPLIER determines the radius of nearby cubes to check.
	set_t tempCubes = getCubes(pt, LOCAL_CUBE_CHECK_MULTIPLIER*res);
	
	///reduce this set to the ones that are highlighted.
	set_t nearbyCubes = alloc_set(0);
	for(i = 0; i<size_set(tempCubes); i++){
		LatticeObj * theObj = objs[tempCubes[i]];
		if( theObj->highlight == HIGHLIGHT_INTERIOR || theObj->highlight == HIGHLIGHT_EXTERIOR ){
			nearbyCubes = put_set(nearbyCubes, tempCubes[i]);
		}
	}
	free_set(tempCubes);

	///clear out some if there are too many	
	while(size_set(nearbyCubes) > LOCAL_CUBE_CHECK_MAX_CUBES){
		int randomVal = random() % size_set(nearbyCubes);
		remove_set(nearbyCubes, nearbyCubes[randomVal]);
	}

	///first we will separate all cubes into HIGHLIGHT_INTERIOR, HIGHLIGHT_EXTERIOR
	set_t insideCubes = alloc_set(SP_MAP);
	set_t outsideCubes = alloc_set(SP_MAP);
	for(i = 0; i<size_set(nearbyCubes); i++){
		//these are the nearby lattice objs.
		LatticeObj * theObj = objs[nearbyCubes[i]];
		if( theObj->highlight == HIGHLIGHT_INTERIOR ){
			insideCubes = associate_set(insideCubes, size_set(insideCubes), theObj);
			continue;
		}
		if( theObj->highlight == HIGHLIGHT_EXTERIOR ){
			outsideCubes = associate_set(outsideCubes, size_set(outsideCubes), theObj);
			continue;
		}
	}
	
	//printf("NUM inside: %i    NUM outside: %i\n", size_set(insideCubes), size_set(outsideCubes) );

	///Make sure that at least one of these sets has nonzero number of cubes, otherwise exit.
	if( (size_set(insideCubes) == 0) && (size_set(outsideCubes) == 0) ){
		//if both lists are empty, return failure.
		//printf("WARNING: LOCAL_CUBE_TEST: fail.\n");
		free_set(insideCubes);
		free_set(outsideCubes);
		free_set(nearbyCubes);
	//	printf("LocalCubes failure: No local flagged cubes to test.\n");
		return 0;
	}

	//Second, we will get all the corner points, nonredundantly, except coaxial points.
	SmallPointHash * insideHash = new SmallPointHash();
	SmallPointHash * outsideHash = new SmallPointHash();

	//list the insideCubes centroids, nonredundantly
	for(i = 0; i<size_set(insideCubes); i++){
	      LatticeObj * thisObj = (LatticeObj*) mapsto_set(insideCubes, i);
	      ///the corners of this cube. (8*3 doubles, for queries and insertion)
	      double * cubePts = thisObj->corners;
	      ///I am going to use the centroid, instead of the corners.
	      double xcent=0, ycent=0, zcent=0;
	      for(j = 0; j<8; j++){
	              xcent += cubePts[3*j+0];
	              ycent += cubePts[3*j+1];
	              zcent += cubePts[3*j+2];
	      }
	      xcent /= 8;     ycent /= 8;     zcent /= 8;
	
	      ///we eliminate anything that is co-axial with pt, for accuracy reasons (low probability).
	      if(pt[0]==xcent || pt[1]==ycent || pt[2]==zcent){continue;}
	              
	      insideHash->addPoint(xcent, ycent, zcent);
	}
	
	//list the outsideCubes centroids, nonredundantly
	for(i = 0; i<size_set(outsideCubes); i++){
	      LatticeObj * thisObj = (LatticeObj*) mapsto_set(outsideCubes, i);
	      ///the centroid of this cube. (8*3 doubles, for queries and insertion)
	      double * cubePts = thisObj->corners;
	      ///I am going to use the centroid, instead of the corners.
	      double xcent=0, ycent=0, zcent=0;
	      for(j = 0; j<8; j++){
	              xcent += cubePts[3*j+0];
	              ycent += cubePts[3*j+1];
	              zcent += cubePts[3*j+2];
	      }
	      xcent /= 8;     ycent /= 8;     zcent /= 8;
	              
	      ///we eliminate anything that is co-axial with pt, for accuracy reasons (low probability).
	      if(pt[0]==xcent || pt[1]==ycent || pt[2]==zcent){continue;}
	              
	      outsideHash->addPoint(xcent, ycent, zcent);
	}
	
	//printf("Total # Points to Test: inner: %i, outer: %i, total: %i\n", insideHash->size(), outsideHash->size(), insideHash->size()+outsideHash->size() );

	///Third, we will begin segment testing.  +1 in the case of zero size.
	int * insideCounts = new int[insideHash->size()+1];
	int * outsideCounts = new int[outsideHash->size()+1];

	///get the intersection count for all inside points
	for(i = 0; i<insideHash->size(); i++){
		//hashPt is passed by reference, do not delete
		double * hashPt = insideHash->getPoint(i);	
		//the final NULL is to ignore distance values, pass a pointer for an array of distance vals
		insideCounts[i] = countIntersectionPts(pt, hashPt, NULL);
	}
	
	///get the intersection count for all outside points
	for(i = 0; i<outsideHash->size(); i++){
		//hashPt is passed by reference, do not delete
		double * hashPt = outsideHash->getPoint(i);	
		//the final NULL is to ignore distance values, pass a pointer for an array of distance vals
		outsideCounts[i] = countIntersectionPts(pt, hashPt, NULL);
	}
	
	///now we look at intersection counts and make a decision.
	//first, check for in-unit contradictions.
	int insideCheck = parityTest(insideCounts, insideHash->size());
	int outsideCheck = parityTest(outsideCounts, outsideHash->size());

	result = -2;	///-2 means we have no answer yet.
	
	//now we check all cases.

	//consistency check with previous (handles case (-1,-1)
	if( insideCheck == -1 && outsideCheck == -1 ){
		printf("WARNING: LOCAL_CUBE_TEST: fail. THIS SHOULD NOT HAPPEN HERE\n");
		printf("WARNING: SHOULD HAVE BEEN CAUGHT EARLIER\n");
		result = 0;
	}

	///fail if either is inconsistent (handles cases (0,-1)(0,0)(0,1)(0,2)(-1,0)(1,0)(2,0) )
	if( insideCheck == 0 || outsideCheck == 0 ){ result = 0; }
	
	///fail if they are inconsistent with each other (handles cases (1,1)(2,2) )
	if( (insideCheck == 2 && outsideCheck == 2) || (insideCheck == 1 && outsideCheck == 1) ){result = 0;}
	
	/// if we want corroborating info from both inside and outside squares, set this to false;
	///This func returns -1 if outside, +1 if inside, 0 if failure.
	bool acceptOneSidedInfo = true;
	///This handles (-1,1)(-1,2)(1,-1)(2,-1)
	if( acceptOneSidedInfo ){
		if( insideCheck == -1 && outsideCheck == 1 ){ result = 1; }		//this is inside
		if( insideCheck == -1 && outsideCheck == 2 ){ result = -1; }	//this is outside

		if( insideCheck == 1 && outsideCheck == -1 ){ result = -1; }	//this is outside
		if( insideCheck == 2 && outsideCheck == -1 ){ result = 1; }		//this is inside
	}
	else{
		if( insideCheck == -1 || outsideCheck == -1 ){ result = 0; }	//if no one sided, fail.
	}

	//handle cases (1,2)(2,1)	
	if( insideCheck == 1 && outsideCheck == 2 ){ result = -1; }
	if( insideCheck == 2 && outsideCheck == 1 ){ result = 1; }


	///clean up
	free_set(insideCubes);
	free_set(outsideCubes);
	free_set(nearbyCubes);
	delete(insideHash);
	delete(outsideHash);
	delete[](insideCounts);
	delete[](outsideCounts);
	
	if(result == -2){
		printf("ERROR: THIS SHOULD NEVER HAPPEN! SOME CASE IN isInsideLocalCubes() was not caught!\n");
		exit(1);
	}
	
	if(result == 0){printf("LOCALCUBES inconsistent (inside: %i, outside: %i)\n", insideCheck, outsideCheck);}
	
	return result;

}
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///This function tests inside/outside based on segment testing to nearby cubes.
///The only cubes tested are those that are highlighted as inside or outside.
///This func returns -1 if outside, +1 if inside, 0 if failure.
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
///################################################################################################///
///# ####   ###   #### ### ####   ## ### #    ##    ##   ### ## ###################################///
///# ### ### # ### ## # ## ### ### # ### # ### # #### ### # #### ##################################///
///# ### ### # ##### ### # ### ##### ### #    ##   ###  ### #### ##################################///
///# ### ### # #####     # ### ##### ### # ### # #######  # #### ##################################///
///# ### ### # ### # ### # ### ### # ### # ### # #### ### # #### ##################################///
///#   ##   ###   ## ### #   ##   ###   ##    ##    ##   ### ## ###################################///
///################################################################################################///
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////








/*
///this determines if the given point is inside the surface
bool TriangleLatticeHash::isInside(double * pt, double * dist)
{
//	////DEBUG
//	bool result;
//	bool rayTest = isInsideRayTest(pt, dist);
//	if(rayTest){result = true;}
//	else{result = false;}
//	return result;
//	////DEBUG


	int i = 0;
	double * tempTri = new double[9];
	for(i = 0; i<3; i++){   
		tempTri[3*i+0] = pt[0];    
		tempTri[3*i+1] = pt[1];    
		tempTri[3*i+2] = pt[2];  
	}
	set_t cubes = getCubes( tempTri );
	delete[](tempTri);

	///if getCubes comes back with size zero, then we are definitely outside.
	if(size_set(cubes) == 0){
		free_set(cubes);
		dist[0] = HUGE_VAL;
		return false;
	}

	///if getCubes is non-empty, we are inside the hash table.
	///find the set of closest triangles
	double * distance = new double[1];
	set_t closeTris = distTo( pt, distance);
	
	/////status provides the current status of analyzing the triangles:
	//-2: no triangles have been analyzed (this stops being the case after first iteration)
	//-1: normals have all been negative so far
	// 0: normals have all been zero so far (this breaks out immidiately)
	// 1: normals have all been positive so far
	//-3: normals are conflicting.  break out and analyze with ray tracing.
	int status = -2;
	
	printf("Number of close Triangles: %i  Distance: %f\n", size_set(closeTris), distance[0]);

	///the current modification is to pass to ray-casting if there is more than one triangle
	///I appear to have errors in multi-triangle calculations that are best left out.
	if(size_set(closeTris) == 1){
	
		for(i = 0; i<size_set(closeTris); i++){
			double * triNorm = new double[3];
			triNorm[0] = surf->triangleNormals[3*closeTris[i]+0];
			triNorm[1] = surf->triangleNormals[3*closeTris[i]+1];
			triNorm[2] = surf->triangleNormals[3*closeTris[i]+2];
		
			double * triPts = new double[3];
			triPts[0] = surf->surfacePoints[3* (surf->triangles[3*closeTris[i]+0]) +0];
			triPts[1] = surf->surfacePoints[3* (surf->triangles[3*closeTris[i]+0]) +1];
			triPts[2] = surf->surfacePoints[3* (surf->triangles[3*closeTris[i]+0]) +2];
		
			double * vector = new double[3];
			vector[0] = triPts[0] - pt[0];
			vector[1] = triPts[1] - pt[1];
			vector[2] = triPts[2] - pt[2];
			
			double test = DOT(triNorm, vector );
			delete[](triNorm); delete[](triPts); delete[](vector);
			
			//printf("Dot product: %f\n", test);
			
			///If the poit is not any surface, throw this error:
			if(test == 0 && distance != 0){ printf("ERROR: NON-MANIFOLD SURFACE GEOMETRY DETECTED AMONG TRIANGLES\n"); exit(1); }
			if(status == -2){
				if(test < 0){ status = -1; }
				if(test > 0){ status = 1; }
				if(test == 0){ status = 0; break;}///here clearly the distance is already zero.
			}
			else{
				if(test == 0){ break; }///here clearly the distance is already zero.
				if(status == -1 && test < 0){ continue; }
				if(status == -1 && test > 0){ status = -3; break; }
				if(status ==  1 && test > 0){ continue; }
				if(status ==  1 && test < 0){ status = -3; break; }
			}		
		}

	}
	else{
		status = -3;
	}

	//if(status == -2){printf("WTF ERROR STATUS IS -2 THIS SHOULD NOT HAPPEN\n");}

	bool result;
	bool rayTest;
	switch(status){
		case -3:
			//printf("beginning ray test\n");
			rayTest = isInsideRayTest(pt, distance);
			if(rayTest){result = true;}
			else{result = false;}
			break;
		case 0:
			result = true;
			break;
		case -1:
			result = false;
			break;
		case 1:
			result = true;
			break;
	}

	if( dist != NULL){ dist[0] = distance[0]; }

	delete[](distance);
	free_set(closeTris);
	free_set(cubes);

	return result;
}
*/




//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///This WRAPPER starts the flood fill
///using flood fill, marks all cubes not containing geometry (Triangles) by setting their highlights to true;
void TriangleLatticeHash::floodFillStart(double x, double y, double z, int highlightCode)
{
	double * coord = new double[3];
	coord[0] = x;
	coord[1] = y;
	coord[2] = z;
	
	floodFillStart( coord, highlightCode );
	delete[](coord);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///This starts the flood fill
///using flood fill, marks all cubes not containing geometry (Triangles) by setting their highlights to true;
void TriangleLatticeHash::floodFillStart(double * coord, int highlightCode )
{
	int i = 0;
//	int j = 0;
	double * tempTri = new double[9];
	for(i = 0; i<3; i++){
		tempTri[3*i+0] = coord[0];
		tempTri[3*i+1] = coord[1];
		tempTri[3*i+2] = coord[2];
	}

	///get the cubes which this point intersects 
	///(it can be more than one if it is on the boundary)
	set_t tempList = getCubes(tempTri);
//	set_t temphash = alloc_set(0);
	delete[](tempTri);
	
	///if this coord intersects with non-empty cubes, or if they are already highlighted, quit.
	for(i = 0; i<size_set(tempList); i++){
		if( (size_set( objs[tempList[i]]->indices ) != 0) || (objs[tempList[i]]->highlight  != NOT_HIGHLIGHTED) ){
			free_set(tempList);
			return;
		}
	}

/*	///list all intersecting cubes that are not empty
	for(i = 0; i<size_set(tempList); i++){
		if( size_set( objs[tempList[i]]->indices ) != 0 ){ temphash = put_set(temphash, tempList[i]); }
	}

	///remove the list of non-empty cubes from the intersecting cubes
	for(i = 0; i<size_set(temphash); i++){
		remove_set(tempList, temphash[i]);
	}
	free_set(temphash);
*/

	///run flood fill on the empty/intersecting cubes remaining 
	floodFillMarker( tempList, highlightCode );

	free_set(tempList);	
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///This continues the floodfill, assuming all latticeIndices are empty.
void TriangleLatticeHash::floodFillMarker( set_t latticeIndices, int highlightCode  )
{
	int i = 0;	int j = 0;

	///create the wavefronts
	set_t wavefront = alloc_set(0);
	for(i = 0; i<size_set(latticeIndices); i++){
		wavefront = put_set(wavefront, latticeIndices[i]);
	}
	set_t nextwave = alloc_set(0);

	int counter;
	int totalCounter = 0;
	///while there is still a frontier
	while( size_set(wavefront) > 0 ){
		counter = 0;

		///For each member of the wavefront, 
		for(i = 0; i<size_set(wavefront); i++){
			int parentObj = wavefront[i];

			///highlight it
			objs[parentObj]->highlight = highlightCode;

			///if its been added to the next wave, remove it from the next wave.
			if(contains_set(nextwave, parentObj)){ remove_set(nextwave, parentObj); }

			//Go through its adjcency list.
			for(j = 0; j<( size_set(objs[parentObj]->adj) ); j++){
				int childObj = objs[parentObj]->adj[j];

				///if the member of the adj list is empty and it is not highlighted,
				///add it to the next wave
				if( (size_set(objs[childObj]->indices) == 0) && (objs[childObj]->highlight == NOT_HIGHLIGHTED) ){
					nextwave = put_set(nextwave, childObj);
					counter++;
				}
			}
		}
		if(counter != 0){
			//printf("Added %i new cubes (out of %i total cubes) in this wave, labelled %i \n", counter, xdim*ydim*zdim, highlightCode );
			totalCounter += counter;
		}

		free_set(wavefront);
		wavefront = nextwave;
		nextwave = alloc_set(0);
	}

	free_set(wavefront);
	free_set(nextwave);

}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///If completely within an exterior-highlighted or interior-highlighted region, 
///this function returns the type of region.  Otherwise it returns NOT_HIGHLIGHTED
int TriangleLatticeHash::insideHighlightedCube(double * coord)
{
	num_highlit_cube_tests++;
	
	int i = 0;
	double * tempTri = new double[9];
	for(i = 0; i<3; i++){
		tempTri[3*i+0] = coord[0];
		tempTri[3*i+1] = coord[1];
		tempTri[3*i+2] = coord[2];
	}

	///get the cubes which this triangle intersects
	set_t tempList = getCubes(tempTri);

	bool allInterior = true;
	bool allExterior = true;
	for( i = 0; i<size_set(tempList); i++){
		if( objs[tempList[i]]->highlight != HIGHLIGHT_INTERIOR ){
			allInterior = false;
			break;
		}
	}
	for( i = 0; i<size_set(tempList); i++){
		if( objs[tempList[i]]->highlight != HIGHLIGHT_EXTERIOR ){
			allExterior = false;
			break;
		}
	}

	delete[](tempTri);
	free_set(tempList);

	int result = NOT_HIGHLIGHTED;
	if(allInterior==true && allExterior==false){ result = HIGHLIGHT_INTERIOR; }
	if(allInterior==false && allExterior==true){ result = HIGHLIGHT_EXTERIOR; }

	return result;
}





//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///eliminates all cube highlights (sets to false)
void TriangleLatticeHash::clearAllHighlights()
{
	int i = 0;
	int j = 0;
	int k = 0;
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				LatticeObj * tempObj = objs[getCube(i, j, k)];

				tempObj->highlight = NOT_HIGHLIGHTED;
			}
		}
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///find Centroid of intersection points between input segment and the surface
double * TriangleLatticeHash::getAverageIntersectionPt(double * a, double * b)
{
	int i = 0;	int j = 0;	int k = 0;	//int l = 0;

	///find the cubes that this segment intersects with.
	double sxmax = a[0]; if(b[0] > a[0]){ sxmax = b[0]; }
	double sxmin = a[0]; if(b[0] < a[0]){ sxmin = b[0]; }
	double symax = a[1]; if(b[1] > a[1]){ symax = b[1]; }
	double symin = a[1]; if(b[1] < a[1]){ symin = b[1]; }
	double szmax = a[2]; if(b[2] > a[2]){ szmax = b[2]; }
	double szmin = a[2]; if(b[2] < a[2]){ szmin = b[2]; }

	int ihi, ilo, jhi, jlo, khi, klo;
	ihi = (int)( (sxmax-xmin)/res );
	jhi = (int)( (symax-ymin)/res );
	khi = (int)( (szmax-zmin)/res );
	if( fmod( (sxmin-xmin), res ) == 0 ){ ilo = ((int)((sxmin-xmin)/res)) -1; } else{ ilo = ((int)((sxmin-xmin)/res)); }
	if( fmod( (symin-ymin), res ) == 0 ){ jlo = ((int)((symin-ymin)/res)) -1; } else{ jlo = ((int)((symin-ymin)/res)); }
	if( fmod( (szmin-zmin), res ) == 0 ){ klo = ((int)((szmin-zmin)/res)) -1; } else{ klo = ((int)((szmin-zmin)/res)); }

	///bound the indices;  These are inclusive indices;
	if( ihi >= xdim ){ ihi = xdim-1; }	if( ilo < 0 ){ ilo = 0; }
	if( jhi >= ydim ){ jhi = ydim-1; }	if( jlo < 0 ){ jlo = 0; }
	if( khi >= zdim ){ khi = zdim-1; }	if( klo < 0 ){ klo = 0; }

//	printf("ilo: %i ihi: %i jlo: %i jhi: %i klo: %i khi: %i  RES: %f\n", ilo, ihi, jlo, jhi, klo, khi, res);
//	printf("ilo: %f ihi: %f jlo: %f jhi: %f klo: %f khi: %f\n", xmin+ilo*res, xmin+ihi*res, ymin+jlo*res, ymin+jhi*res, zmin+klo*res, zmin+khi*res);
//	printf("length: %f\n", vectorSize(b[0]-a[0], b[1]-a[1], b[2]-a[2]));

	///get a list of the cubes within the bounds
	set_t tempCubes = alloc_set(0);
	for(i = ilo; i<=ihi; i++){
		for(j = jlo; j<=jhi; j++){
			for(k = klo; k<=khi; k++){
				LatticeObj * tempObj = objs[ grid[i][j][k] ];

				///test to see if these cubes actually intersect the segment.
				if( tempObj->segmentIntersect( a,b ) ){
					tempCubes = put_set(tempCubes, grid[i][j][k] );
				}
			}
		}
	}

	//printf("number of cubes inspected for segment/triangle intersection: %i\n", size_set(tempCubes));

	///accumulate the intersection points
	double * resultPts = new double[1];  resultPts[0] = 0;
	for(i = 0; i<size_set(tempCubes); i++){
		LatticeObj * tempObj = objs[ tempCubes[i] ];
		double * tempResult = tempObj->getIntersectPts( a, b );
		
		///print out cube coords:
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*0+0], tempObj->corners[3*0+1], tempObj->corners[3*0+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*1+0], tempObj->corners[3*1+1], tempObj->corners[3*1+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*2+0], tempObj->corners[3*2+1], tempObj->corners[3*2+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*3+0], tempObj->corners[3*3+1], tempObj->corners[3*3+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*4+0], tempObj->corners[3*4+1], tempObj->corners[3*4+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*5+0], tempObj->corners[3*5+1], tempObj->corners[3*5+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*6+0], tempObj->corners[3*6+1], tempObj->corners[3*6+2] );
	//	printf("Cube: %f %f %f\n", tempObj->corners[3*7+0], tempObj->corners[3*7+1], tempObj->corners[3*7+2] );

		//printf("Number of intersectionPts Found in this cube: %i  Num triangles Tested: %i\n", (int)tempResult[0], size_set(tempObj->indices));
		if(tempResult[0] > 0){
			double * newResult = new double[1 + (3*((int)resultPts[0])) + (3*((int)tempResult[0]))];
			for(j = 0; j<resultPts[0]; j++){
				newResult[1+(3*j+0)] = resultPts[1+(3*j+0)];
				newResult[1+(3*j+1)] = resultPts[1+(3*j+1)];
				newResult[1+(3*j+2)] = resultPts[1+(3*j+2)];
			}
			for(j = 0; j<tempResult[0]; j++){
				newResult[1+(3*((int)resultPts[0]))+(3*j+0)] = tempResult[1+(3*j+0)];
				newResult[1+(3*((int)resultPts[0]))+(3*j+1)] = tempResult[1+(3*j+1)];
				newResult[1+(3*((int)resultPts[0]))+(3*j+2)] = tempResult[1+(3*j+2)];
			}
			newResult[0] = resultPts[0] + tempResult[0];
			delete[](resultPts);
			resultPts = newResult;
		}
		delete[](tempResult); ///this is never null, so always delete.
	}
	
	///average the intersection points
	double * result = NULL;
	if(resultPts[0]==0){
#ifdef FORCE_LOGICAL_CONSISTANCY
		double d1,d2,t;
		set_t extraSet;
		extraSet = distTo(a, &d1);	free_set(extraSet);
		extraSet = distTo(b, &d2);	free_set(extraSet);
		if((d1<res && d2<res) && d1<d2){t = .25; }
		if((d1<res && d2<res) && d2<d1){t = .75; }
		else{t = .5;}

		delete[](resultPts);
		resultPts = new double[4];
		resultPts[0] = 1;
		resultPts[1] = a[0] + t*(b[0]-a[0]);
		resultPts[2] = a[1] + t*(b[1]-a[1]);
		resultPts[3] = a[2] + t*(b[2]-a[2]);
		printf("WARNING: EDGE CONSISTANCY FORCED\n");		
		printf("a[0]=%f; a[1]=%f; a[2]=%f;\n", a[0], a[1], a[2]);
		printf("b[0]=%f; b[1]=%f; b[2]=%f;\n", b[0], b[1], b[2]);
		printf("pt[0]=%f; pt[1]=%f; pt[2]=%f;\n", resultPts[1], resultPts[2], resultPts[3]);
#else 
		delete[](resultPts);
		printf("ERROR: AVERAGE INTERSECTION POINT IN TRIANGLE HASH [%s] IS NULL!\n", name);
		printf("a[0]=%f; a[1]=%f; a[2]=%f;\n", a[0], a[1], a[2]);
		printf("b[0]=%f; b[1]=%f; b[2]=%f;\n", b[0], b[1], b[2]);

		///isInside will call isInsideRayTest if it becomes indecisive.
		bool test1 = isInside( a );
		bool test2 = isInside( b );

		if(test1 && !test2){ printf("ERROR: HOWEVER A IS INSIDE AND B IS OUTSIDE\n"); } // exit(1); }
		if(!test1 && test2){ printf("ERROR: HOWEVER A IS OUTSIDE AND B IS INSIDE\n"); } // exit(1); }
		if(test1 && test2){  printf("ERROR: HOWEVER BOTH POINTS ARE INSIDE\n"); } // exit(1); }
		if(!test1 && !test2){printf("ERROR: HOWEVER BOTH POINTS ARE OUTSIDE\n"); } // exit(1); }
		exit(1);
#endif

	}

	result = new double[3];
	result[0] = 0;	result[1] = 0;	result[2] = 0;
	for(i = 0; i<resultPts[0]; i++){
		result[0] += resultPts[1+(3*i+0)];
		result[1] += resultPts[1+(3*i+1)];
		result[2] += resultPts[1+(3*i+2)];
	}
	result[0] /= resultPts[0];	result[1] /= resultPts[0];	result[2] /= resultPts[0];

//	printf("A: %f %f %f    AvgInter: %f %f%f     B: %f %f %f\n", a[0], a[1], a[2], result[0], result[1], result[2], b[0], b[1], b[2] );



	///output
	delete[](resultPts);
	return result;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///This is the same as above but does not demand that the poitns be inside and outside - so if there is no intersection, it returns NULL.
///find Centroid of intersection points between input segment and the surface
double * TriangleLatticeHash::getAverageIntersectionPtNoTest(double * a, double * b)
{
	int i = 0;	int j = 0;	int k = 0;	//int l = 0;

	///find the cubes that this segment intersects with.
	double sxmax = a[0]; if(b[0] > a[0]){ sxmax = b[0]; }
	double sxmin = a[0]; if(b[0] < a[0]){ sxmin = b[0]; }
	double symax = a[1]; if(b[1] > a[1]){ symax = b[1]; }
	double symin = a[1]; if(b[1] < a[1]){ symin = b[1]; }
	double szmax = a[2]; if(b[2] > a[2]){ szmax = b[2]; }
	double szmin = a[2]; if(b[2] < a[2]){ szmin = b[2]; }

	int ihi, ilo, jhi, jlo, khi, klo;
	ihi = (int)( (sxmax-xmin)/res );
	jhi = (int)( (symax-ymin)/res );
	khi = (int)( (szmax-zmin)/res );
	if( fmod( (sxmin-xmin), res ) == 0 ){ ilo = ((int)((sxmin-xmin)/res)) -1; } else{ ilo = ((int)((sxmin-xmin)/res)); }
	if( fmod( (symin-ymin), res ) == 0 ){ jlo = ((int)((symin-ymin)/res)) -1; } else{ jlo = ((int)((symin-ymin)/res)); }
	if( fmod( (szmin-zmin), res ) == 0 ){ klo = ((int)((szmin-zmin)/res)) -1; } else{ klo = ((int)((szmin-zmin)/res)); }

	///bound the indices;  These are inclusive indices;
	if( ihi >= xdim ){ ihi = xdim-1; }	if( ilo < 0 ){ ilo = 0; }
	if( jhi >= ydim ){ jhi = ydim-1; }	if( jlo < 0 ){ jlo = 0; }
	if( khi >= zdim ){ khi = zdim-1; }	if( klo < 0 ){ klo = 0; }

	///get a list of the cubes within the bounds
	set_t tempCubes = alloc_set(0);
	for(i = ilo; i<=ihi; i++){
		for(j = jlo; j<=jhi; j++){
			for(k = klo; k<=khi; k++){
				LatticeObj * tempObj = objs[ grid[i][j][k] ];

				///test to see if these cubes actually intersect the segment.
				if( tempObj->segmentIntersect( a,b ) ){
					tempCubes = put_set(tempCubes, grid[i][j][k] );
				}
			}
		}
	}

	//printf("number of cubes inspected for segment/triangle intersection: %i\n", size_set(tempCubes));

	///accumulate the intersection points
	double * resultPts = new double[1];  resultPts[0] = 0;
	for(i = 0; i<size_set(tempCubes); i++){
		LatticeObj * tempObj = objs[ tempCubes[i] ];
		double * tempResult = tempObj->getIntersectPts( a, b );

		//printf("Number of intersectionPts Found in this cube: %i  Num triangles Tested: %i\n", (int)tempResult[0], size_set(tempObj->indices));
		if(tempResult[0] > 0){
			double * newResult = new double[1 + (3*((int)resultPts[0])) + (3*((int)tempResult[0]))];
			for(j = 0; j<resultPts[0]; j++){
				newResult[1+(3*j+0)] = resultPts[1+(3*j+0)];
				newResult[1+(3*j+1)] = resultPts[1+(3*j+1)];
				newResult[1+(3*j+2)] = resultPts[1+(3*j+2)];
			}
			for(j = 0; j<tempResult[0]; j++){
				newResult[1+(3*((int)resultPts[0]))+(3*j+0)] = tempResult[1+(3*j+0)];
				newResult[1+(3*((int)resultPts[0]))+(3*j+1)] = tempResult[1+(3*j+1)];
				newResult[1+(3*((int)resultPts[0]))+(3*j+2)] = tempResult[1+(3*j+2)];
			}
			newResult[0] = resultPts[0] + tempResult[0];
			delete[](resultPts);
			resultPts = newResult;
		}
		delete[](tempResult); ///this is never null, so always delete.
	}
	
	///average the intersection points
	double * result = NULL;
	if(resultPts[0]>0){
		result = new double[3];
		result[0] = 0;	result[1] = 0;	result[2] = 0;
		for(i = 0; i<resultPts[0]; i++){
			result[0] += resultPts[1+(3*i+0)];
			result[1] += resultPts[1+(3*i+1)];
			result[2] += resultPts[1+(3*i+2)];
		}
		result[0] /= resultPts[0];	result[1] /= resultPts[0];	result[2] /= resultPts[0];
	}

	///output
	delete[](resultPts);
	return result;
}








//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///counts the number of times a segment intersects with triangles from the surface
///dist[0] provides the distance from a to the closest intersection point
int TriangleLatticeHash::countIntersectionPts(double * a, double * b, double * distanceOutput)
{
	int i = 0;
	int j = 0;
	int k = 0;
//	int l = 0;
	
	double shortestDistance = vectorSize(a,b);

	///find the cubes that this segment intersects with.
	double sxmax = a[0]; if(b[0] > a[0]){ sxmax = b[0]; }
	double sxmin = a[0]; if(b[0] < a[0]){ sxmin = b[0]; }
	double symax = a[1]; if(b[1] > a[1]){ symax = b[1]; }
	double symin = a[1]; if(b[1] < a[1]){ symin = b[1]; }
	double szmax = a[2]; if(b[2] > a[2]){ szmax = b[2]; }
	double szmin = a[2]; if(b[2] < a[2]){ szmin = b[2]; }

	int ihi, ilo, jhi, jlo, khi, klo;
	ihi = (int)( (sxmax-xmin)/res );
	jhi = (int)( (symax-ymin)/res );
	khi = (int)( (szmax-zmin)/res );
	if( fmod( (sxmin-xmin), res ) == 0 ){ ilo = ((int)((sxmin-xmin)/res)) -1; } else{ ilo = ((int)((sxmin-xmin)/res)); }
	if( fmod( (symin-ymin), res ) == 0 ){ jlo = ((int)((symin-ymin)/res)) -1; } else{ jlo = ((int)((symin-ymin)/res)); }
	if( fmod( (szmin-zmin), res ) == 0 ){ klo = ((int)((szmin-zmin)/res)) -1; } else{ klo = ((int)((szmin-zmin)/res)); }

	///bound the indices;  These are inclusive indices;
	if( ihi >= xdim ){ ihi = xdim-1; }	if( ilo < 0 ){ ilo = 0; }
	if( jhi >= ydim ){ jhi = ydim-1; }	if( jlo < 0 ){ jlo = 0; }
	if( khi >= zdim ){ khi = zdim-1; }	if( klo < 0 ){ klo = 0; }

//	printf("ilo: %i ihi: %i jlo: %i jhi: %i klo: %i khi: %i  RES: %f\n", ilo, ihi, jlo, jhi, klo, khi, res);
//	printf("ilo: %f ihi: %f jlo: %f jhi: %f klo: %f khi: %f\n", xmin+ilo*res, xmin+ihi*res, ymin+jlo*res, ymin+jhi*res, zmin+klo*res, zmin+khi*res);
//	printf("length: %f\n", vectorSize(b[0]-a[0], b[1]-a[1], b[2]-a[2]));

	///get a list of the cubes within the bounds
	set_t tempCubes = alloc_set(0);
	for(i = ilo; i<=ihi; i++){
		for(j = jlo; j<=jhi; j++){
			for(k = klo; k<=khi; k++){
				LatticeObj * tempObj = objs[ grid[i][j][k] ];

				///test to see if these cubes actually intersect the segment.
				if( tempObj->segmentIntersect( a,b ) ){
					tempCubes = put_set(tempCubes, grid[i][j][k] );
				}
			}
		}
	}
	
//	printf("number of cubes inspected for segment/triangle intersection: %i\n", size_set(tempCubes));

	///accumulate the intersection points, eliminating duplicates.
	set_t iPoints = alloc_set(SP_MAP);
	for(i = 0; i<size_set(tempCubes); i++){
		LatticeObj * tempObj = objs[ tempCubes[i] ];
		double * tempResult = tempObj->getIntersectPts( a, b );

		if(tempResult[0] > 0){
			double * thisPt = new double[3];
			for(j = 0; j<tempResult[0]; j++){
				thisPt[0] = tempResult[1+(3*j+0)];
				thisPt[1] = tempResult[1+(3*j+1)];
				thisPt[2] = tempResult[1+(3*j+2)];
				
	//			printf("Intersect: %f %f %f tempObj#: %i\n", thisPt[0], thisPt[1], thisPt[2], tempCubes[i] );
				
				bool okForInsertion = true;
				///if there is nothing too similar to these points, insert.
				for(k = 0; k<size_set(iPoints); k++){
					double * oldPt = (double *) mapsto_set(iPoints, k);
					double dist = vectorSize(thisPt, oldPt);
					if(dist < .00001){
						okForInsertion = false;
						break;
					}
				}
				if(okForInsertion){
					double * insert = new double[3];
					insert[0] = thisPt[0]; insert[1] = thisPt[1]; insert[2] = thisPt[2];
					iPoints = associate_set(iPoints, size_set(iPoints), insert);
					double tempDist = vectorSize(insert, a);
					if(tempDist < shortestDistance){
						shortestDistance = tempDist;
					}
				}
			}
			delete[](thisPt);
		}
		delete[](tempResult); ///this is never null, so always delete.
	}
	
	int result = size_set(iPoints);
	
	for(i = 0; i<size_set(iPoints); i++){
		double * insert = (double *) mapsto_set(iPoints, i);
		delete[](insert);
	}
	free_set(iPoints);
	free_set(tempCubes);

	if(distanceOutput != NULL){ distanceOutput[0] = shortestDistance; }

	return result;
}







/////////////////////////////////////////////////////////////////////////////
///Benchmarking function counts the number of flagged points, unflagged points
///  in the entire lattice of corners.
//	a tallied point gets zero if it has no flag, -1 if it's outside, +1 if inside.
/////////////////////////////////////////////////////////////////////////////
///this code uses the following cube layout, with 3D axes to the right.
///      7---6            ^
///      |\  |\       \   |
///      | 4---5       +y |
///      | | | |        \ z+
///      3-|-2 |         \|  
///       \|  \|          0--x+-->    
///        0---1
/////////////////////////////////////////////////////////////////////////////
void TriangleLatticeHash::surveyFlagStatus()
{
	printf("Counting interior/exterior points from the triangle lattice hash...\n");
	printf("xdim: %i\n", xdim);
	printf("ydim: %i\n", ydim);
	printf("zdim: %i\n", zdim);
	
	int i = 0;
	int j = 0;
	int k = 0;
	//int maxPoints = (xdim+1)*(ydim+1)*(zdim+1)

	int *** tally = new int**[xdim+1];

	/// set up the tally array.
	for(i = 0; i<xdim+1; i++){
		tally[i] = new int*[ydim+1];
		for(j = 0; j<ydim+1; j++){
			tally[i][j] = new int[zdim+1];
			for(k = 0; k<zdim+1; k++){
				tally[i][j][k] = 0;
			}
		}
	}
	
	///iterate through the grid array, which indexes the lattice cubes.
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				
				//find out how the object is flagged.
				int gridIndex = grid[i][j][k];
				LatticeObj * tempObj = objs[gridIndex];
				
				//if the object is flagged interior, mark the tally interior.
				if(tempObj->highlight == HIGHLIGHT_INTERIOR){
					tally[i+0][j+0][k+0] = 1;
					tally[i+0][j+0][k+1] = 1;
					tally[i+0][j+1][k+0] = 1;
					tally[i+0][j+1][k+1] = 1;
					tally[i+1][j+0][k+0] = 1;
					tally[i+1][j+0][k+1] = 1;
					tally[i+1][j+1][k+0] = 1;
					tally[i+1][j+1][k+1] = 1;
				}
			
				//if the object is flagged exterior, mark the tally exterior.
				if(tempObj->highlight == HIGHLIGHT_EXTERIOR){
					tally[i+0][j+0][k+0] = -1;
					tally[i+0][j+0][k+1] = -1;
					tally[i+0][j+1][k+0] = -1;
					tally[i+0][j+1][k+1] = -1;
					tally[i+1][j+0][k+0] = -1;
					tally[i+1][j+0][k+1] = -1;
					tally[i+1][j+1][k+0] = -1;
					tally[i+1][j+1][k+1] = -1;
				}
			}
		}
	}
	
	///count the tally.
	int numInterior = 0;
	int numExterior = 0;
	int numUnmarked = 0;
	
	for(i = 0; i<xdim+1; i++){
		for(j = 0; j<ydim+1; j++){
			for(k = 0; k<zdim+1; k++){
				int myVal = tally[i][j][k];
				if( myVal == -1){ numExterior++; }
				if( myVal ==  1){ numInterior++; }
				if( myVal ==  0){ numUnmarked++; }
			}
		}
	}

	printf("Tally data from Triangle Lattice Hash:\n");
	printf("Number of Interior Points: %i\n", numInterior);
	printf("Number of Exterior Points: %i\n", numExterior);
	printf("Number of Unmarked Points: %i\n", numUnmarked);

	///clear the data
	for(i = 0; i<xdim+1; i++){
		for(j = 0; j<ydim+1; j++){
			delete[](tally[i][j]);
		}
		delete[](tally[i]);
	}
	delete[](tally);
}



/////////////////////////////////////////////////////////////////////////////
///Benchmarking function counts how many calls were made from external query.
void TriangleLatticeHash::surveyTestStatus()
{
	printf("Test data from Triangle Lattice Hash:\n");
	printf("Number of highlit cube tests: %i\n", num_highlit_cube_tests);
	printf("Number of cache hits: %i\n", num_cache_hits );
	printf("Number of local cubes tests: %i\n", num_local_cubes_tests);
	printf("Number of long ray tests: %i\n", num_long_ray_tests);

}




/////////////////////////////////////////////////////////////////////////////
//This function tells you if a given triangle intersects with a triangle in the surface
bool TriangleLatticeHash::triangleIntersection( double * ta, double * tb, double * tc )
{
	int i = 0;
	int j = 0;
	int result = 0;
	
	double * triangle = new double[9];
	triangle[0] = ta[0];
	triangle[1] = ta[1];
	triangle[2] = ta[2];
	triangle[3] = tb[0];
	triangle[4] = tb[1];
	triangle[5] = tb[2];
	triangle[6] = tc[0];
	triangle[7] = tc[1];
	triangle[8] = tc[2];

	double * tri0 = new double[3];
	double * tri1 = new double[3];
	double * tri2 = new double[3];

	///returns the cubes which this triangle intersects
	set_t myCubes = getCubes(triangle);
	set_t myTriangles = alloc_set(0);
	for(i = 0; i<size_set(myCubes); i++){
		LatticeObj * tempObj = objs[myCubes[i]];
		if( tempObj->triangleIntersect( triangle ) ){
			for(j = 0; j<size_set(tempObj->indices); j++){
				myTriangles = put_set(myTriangles, tempObj->indices[j]);
			}
		}
	}
	free_set(myCubes);

	for(i = 0; i<size_set(myTriangles); i++){
		int myTriIndex = myTriangles[i];
		double * testTriangle = surf->getTriangle(myTriIndex); 
		tri0[0] = testTriangle[0];
		tri0[1] = testTriangle[1];
		tri0[2] = testTriangle[2];
		tri1[0] = testTriangle[3];
		tri1[1] = testTriangle[4];
		tri1[2] = testTriangle[5];
		tri2[0] = testTriangle[6];
		tri2[1] = testTriangle[7];
		tri2[2] = testTriangle[8];
		delete[](testTriangle);

		result = tri_tri_overlap_test_3d(ta, tb, tc, tri0, tri1, tri2);
		
		if(result==1){
			break;
		}
		
	}
	
	delete[](triangle);
	free_set(myTriangles);
	delete[](tri0);
	delete[](tri1);
	delete[](tri2);
	
	return result;
}


























