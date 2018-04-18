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
 * File: SurfaceObject.cpp
 *       Implementation of a container for surface data
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

//////////////////////////////////////////////////////////////////////
// SurfaceObject.cpp: implementation of the SurfaceObject class.
//
//  This a class that represents a connolly surface
//  using vectors, norms, and a triangle topology.
//
//  Implementation Note:
//  The version provided below is an abbreviation of the more complex class
//  called SurfaceObject in SurfaceExtractor, which handles colored surfaces
//  and several other functionalities.  Here, this datastructure is used only
//  as a container.
//////////////////////////////////////////////////////////////////////

#include "SurfaceObject.h"


//////////////////////////////////////////////////////////////////////
// Construction
//
// This is the null constructor, which generates no surface.
SurfaceObject::SurfaceObject()
{
        numPoints = 0;
        numTriangles = 0;
        highlights = NULL;
        edges = NULL;
        colors = NULL;

        centroid = new double[1];
        surfacePoints = new double[1];
        surfaceNormals = new double[1];
        triangles = new int[1];
        triangleNormals = new double[1];
}


//////////////////////////////////////////////////////////////////////
// Construction
//
// Construct the object with just a set_t of double[3] points, and 
// a set_t of int[3] triangles.  Indices in the int[3] of triangles
// reference array positions in the set of point vectors.
//////////////////////////////////////////////////////////////////////
SurfaceObject::SurfaceObject(set_t pts, set_t tris)
{
	int i  = 0;
	int j  = 0;
	
	///get the sizes
	numPoints = size_set(pts);
	numTriangles = size_set(tris);
	highlights = NULL;
	edges = NULL;
	colors = NULL;

	///store which triangles the points are attached to, for poitn normals later.
	set_t adjTriangles = alloc_set(SP_MAP);
	for(i = 0; i<numPoints; i++){
		set_t temp = alloc_set(0);
		adjTriangles = associate_set(adjTriangles, i, temp);
	}

	/////fill in the pts (+1 to avoid 0 allocation when an empty surface comes in)
	surfacePoints = new double[3*numPoints+1];
	for(i = 0; i<numPoints; i++){
		//get the pt and store it.
		double * pt = (double *) mapsto_set(pts, pts[i]);
		surfacePoints[3*i+0] = pt[0];
		surfacePoints[3*i+1] = pt[1];
		surfacePoints[3*i+2] = pt[2];
	}

	///fill in the triangles (+1 to avoid 0 allocation with empty surfs)
	triangles = new int[3*numTriangles+1];
	triangleNormals = new double[3*numTriangles+1];
	for(i = 0; i<numTriangles; i++){
		//get the triangle an store it
		int * t = (int *) mapsto_set(tris, tris[i]);
		triangles[3*i+0] = t[0];
		triangles[3*i+1] = t[1];
		triangles[3*i+2] = t[2];
		
		///make this triangle adjacent to it's points.
		set_t myPt;
		myPt = (set_t) mapsto_set(adjTriangles, t[0]);
		//printf("T[0]: %i\n", t[0]);
		myPt = put_set(myPt, i);
		adjTriangles = associate_set(adjTriangles, t[0], myPt);
		myPt = (set_t) mapsto_set(adjTriangles, t[1]);
		myPt = put_set(myPt, i);
		adjTriangles = associate_set(adjTriangles, t[1], myPt);
		myPt = (set_t) mapsto_set(adjTriangles, t[2]);
		myPt = put_set(myPt, i);
		adjTriangles = associate_set(adjTriangles, t[2], myPt);
		
		//now compute normals.
		//get point Vals (do not delete these)
		double * firstPt = (double *) mapsto_set(pts, t[0]);
		double * secndPt = (double *) mapsto_set(pts, t[1]);
		double * thirdPt = (double *) mapsto_set(pts, t[2]);
		
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

		//store triangle normals
		triangleNormals[3*i+0] = normalizedNormal[0];
		triangleNormals[3*i+1] = normalizedNormal[1];
		triangleNormals[3*i+2] = normalizedNormal[2];
		
		////clear and increment
		delete[](first);
		delete[](second);
		delete[](normal);
		delete[](normalizedNormal);
	}
	
	//now that the triangle normals are filled in, now do the point normals.
	surfaceNormals = new double[3*numPoints+1];
	for(i = 0; i<numPoints; i++){
		//zero out the vector
		surfaceNormals[3*i+0] = 0;	surfaceNormals[3*i+1] = 0;	surfaceNormals[3*i+2] = 0;
		set_t myAdj = (set_t) mapsto_set(adjTriangles, i);
		//average the normals
		for(j = 0; j<size_set(myAdj); j++){
			surfaceNormals[3*i+0] += triangleNormals[3*myAdj[j]+0];
			surfaceNormals[3*i+1] += triangleNormals[3*myAdj[j]+1];
			surfaceNormals[3*i+2] += triangleNormals[3*myAdj[j]+2];
		}
		surfaceNormals[3*i+0] /= (double) size_set(myAdj);
		surfaceNormals[3*i+1] /= (double) size_set(myAdj);
		surfaceNormals[3*i+2] /= (double) size_set(myAdj);
	}

	for(i = 0; i<numPoints; i++){
		set_t temp = (set_t) mapsto_set(adjTriangles, i);
		free_set(temp);
	}
	free_set(adjTriangles);

	///compute the centroid
	centroid = getCentroid();

}


//////////////////////////////////////////////////////////////////////
// Construction
//
// This is the standard constructor, designed for file parsing.
// Points and normals are parsed in the same large array, to facilitate parsing.
//////////////////////////////////////////////////////////////////////
SurfaceObject::SurfaceObject(int numPts, double * ptsAndNorms, int numTris, int * inputTris)
{
	///get the sizes
	numPoints = numPts;
	numTriangles = numTris;
	highlights = NULL;
	edges = NULL;
	colors = NULL;
	
	int * tris;

	tris = inputTris;

	///allocate the data (+1 to avoid 0 allocation when an empty surface comes in)
	surfacePoints = new double[3*numPoints+1];
	surfaceNormals = new double[3*numPoints+1];
	triangles = new int[3*numTriangles+1];

	int i = 0;
	///map in the data
	for(i = 0; i<numPoints; i++){
		surfacePoints[3*i+0] = ptsAndNorms[6*i+0];
		surfacePoints[3*i+1] = ptsAndNorms[6*i+1];
		surfacePoints[3*i+2] = ptsAndNorms[6*i+2];
		surfaceNormals[3*i+0] = ptsAndNorms[6*i+3];
		surfaceNormals[3*i+1] = ptsAndNorms[6*i+4];
		surfaceNormals[3*i+2] = ptsAndNorms[6*i+5];
	}

	////compute the centroid
	centroid = getCentroid();
	
	///map in the triangles
	for(i = 0; i<numTriangles; i++){
		triangles[3*i+0] = tris[3*i+0];
		triangles[3*i+1] = tris[3*i+1];
		triangles[3*i+2] = tris[3*i+2];
	}

	///Compute the triangleNormals (+1 to avoid 0 allocation when an empty surface comes in)
	triangleNormals = new double[3*numTriangles+1];
	double * cross = new double[3];
	double * p0 = new double[3];
	double * p1 = new double[3];
	double * p2 = new double[3];
	double * v1 = new double[3];
	double * v2 = new double[3];
	double * tempNormal = new double[3];
	for(i = 0; i<numTriangles; i++){
		p0[0] = surfacePoints[3*triangles[3*i+0]+0];
		p0[1] = surfacePoints[3*triangles[3*i+0]+1];
		p0[2] = surfacePoints[3*triangles[3*i+0]+2];

		p1[0] = surfacePoints[3*triangles[3*i+1]+0];
		p1[1] = surfacePoints[3*triangles[3*i+1]+1];
		p1[2] = surfacePoints[3*triangles[3*i+1]+2];

		p2[0] = surfacePoints[3*triangles[3*i+2]+0];
		p2[1] = surfacePoints[3*triangles[3*i+2]+1];
		p2[2] = surfacePoints[3*triangles[3*i+2]+2];

		v1[0] = p1[0]-p0[0];	v1[1] = p1[1]-p0[1];	v1[2] = p1[2] - p0[2];
		v2[0] = p2[0]-p0[0];	v2[1] = p2[1]-p0[1];	v2[2] = p2[2] - p0[2];
		
		tempNormal[0] = (surfaceNormals[3*triangles[3*i+0]+0] + surfaceNormals[3*triangles[3*i+1]+0] + surfaceNormals[3*triangles[3*i+2]+0])/3;
		tempNormal[1] = (surfaceNormals[3*triangles[3*i+0]+1] + surfaceNormals[3*triangles[3*i+1]+1] + surfaceNormals[3*triangles[3*i+2]+1])/3;
		tempNormal[2] = (surfaceNormals[3*triangles[3*i+0]+2] + surfaceNormals[3*triangles[3*i+1]+2] + surfaceNormals[3*triangles[3*i+2]+2])/3;

		CROSS(cross,v1,v2);
		if( DOT(cross, tempNormal)<0 ){	///if the vectors face opposite way, reverse the x-product
			
		//	printf("WARNING: (Triangle %i) Normals deviate from averaged normals (%i %i %i)!\n", i, triangles[3*i+0], triangles[3*i+1], triangles[3*i+2]);
		//	printf("TEMPNORMAL: %lf %lf %lf (point: %i)\n", tempNormal[0], tempNormal[1], tempNormal[2], triangles[3*i+0] );
		//	printf("CROSS: %lf %lf %lf\n", cross[0], cross[1], cross[2]);
			
			CROSS(cross,v2,v1);
			//printf("CROSS: %lf %lf %lf\n", cross[0], cross[1], cross[2]);
		}
		
		double * tempVec = normalizeVector(cross);
		
		triangleNormals[3*i+0] = tempVec[0];
		triangleNormals[3*i+1] = tempVec[1];
		triangleNormals[3*i+2] = tempVec[2];
		delete[](tempVec);
	}
	delete[](cross); delete[](tempNormal);
	delete[](p0); 	delete[](p1); 	delete[](p2); 
	delete[](v1); 	delete[](v2); 
}


//////////////////////////////////////////////////////////////////////
// Destruction
//////////////////////////////////////////////////////////////////////
SurfaceObject::~SurfaceObject()
{
	if(highlights != NULL){ free_set(highlights); }
	if(edges != NULL){ free_set(edges); }
	if(colors != NULL){ delete[](colors); }

	delete[](surfacePoints);
	delete[](surfaceNormals);
	delete[](triangles);
	delete[](centroid);
	delete[](triangleNormals);
}




//////////////////////////////////////////////////////////////////////
// returns the coords of the requested triangle
// returns null if out of range.
// notice that the triangle must be formatted in this way in the array (3 then 3 then 3)
//////////////////////////////////////////////////////////////////////
double * SurfaceObject::getTriangle(int num)
{
	if( (num < 0) || (num>=numTriangles) ){
		printf("ERROR: Triangle Not Found!\n");
		return NULL;
	}
	
	double * result = new double[9];
	
	result[0] = surfacePoints[ 3*triangles[3*num+0]+0 ];
	result[1] = surfacePoints[ 3*triangles[3*num+0]+1 ];
	result[2] = surfacePoints[ 3*triangles[3*num+0]+2 ];
	result[3] = surfacePoints[ 3*triangles[3*num+1]+0 ];
	result[4] = surfacePoints[ 3*triangles[3*num+1]+1 ];
	result[5] = surfacePoints[ 3*triangles[3*num+1]+2 ];
	result[6] = surfacePoints[ 3*triangles[3*num+2]+0 ];
	result[7] = surfacePoints[ 3*triangles[3*num+2]+1 ];
	result[8] = surfacePoints[ 3*triangles[3*num+2]+2 ];
	
	return result;
}


//////////////////////////////////////////////////////////////////////
// Computes the centroid based on the point positions
//////////////////////////////////////////////////////////////////////
double * SurfaceObject::getCentroid()
{
	int i = 0;
	double * cent = new double[3];
	cent[0] = 0;
	cent[1] = 1;
	cent[2] = 2;
	
	for(i = 0; i<numPoints; i++){
		cent[0] += surfacePoints[3*i+0];
		cent[1] += surfacePoints[3*i+1];
		cent[2] += surfacePoints[3*i+2];
	}

	cent[0] /= numPoints;
	cent[1] /= numPoints;
	cent[2] /= numPoints;

	return cent;
}






//////////////////////////////////////////////////////////////////////
// This function flips the normals of the surfaceObject backwards 
// so that we can treat it as a negative volume
//////////////////////////////////////////////////////////////////////
void SurfaceObject::flipNormals()
{
	int i = 0;
	for(i = 0; i<3*numPoints; i++){
		surfaceNormals[i] = -surfaceNormals[i];
	}
	for(i = 0; i<3*numTriangles; i++){
		triangleNormals[i] = -triangleNormals[i];
	}
	int temp;
	for(i = 0; i<numTriangles; i++){
		triangles[3*i+0] = triangles[3*i+0];
		temp = triangles[3*i+1];
		triangles[3*i+1] = triangles[3*i+2];
		triangles[3*i+2] = temp;		
	}
}






///////////////////////////////////////////////////////
///copies a surfaceObject
SurfaceObject * SurfaceObject::copy()
{
	int i = 0;

	double * ptsAndNorms = new double[6*numPoints];
	for(i = 0; i<numPoints; i++){
		ptsAndNorms[i*6+0] = surfacePoints[3*i+0];
		ptsAndNorms[i*6+1] = surfacePoints[3*i+1];
		ptsAndNorms[i*6+2] = surfacePoints[3*i+2];
		ptsAndNorms[i*6+3] = surfaceNormals[3*i+0];
		ptsAndNorms[i*6+4] = surfaceNormals[3*i+1];
		ptsAndNorms[i*6+5] = surfaceNormals[3*i+2];
	}

	SurfaceObject * result = new SurfaceObject(numPoints, ptsAndNorms, numTriangles, triangles);

	delete[](ptsAndNorms);
	
	return result;
}
///////////////////////////////////////////////////////





/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/




