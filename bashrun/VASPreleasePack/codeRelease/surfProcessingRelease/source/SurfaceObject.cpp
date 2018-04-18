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
 *       Object for representing triangle meshes
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
// This is the standard constructor, designed for file parsing.
// Points and normals are parsed in the same large array, to facilitate parsing.
SurfaceObject::SurfaceObject(int numPts, double * ptsAndNorms, int numTris, int * inputTris, bool elimIntCavities)
{
	int i = 0;

	///get the sizes
	numPoints = numPts;
	numTriangles = numTris;
	highlights = NULL;
	edges = NULL;
	colors = NULL;
	
	int * tris = inputTris;

	///allocate the data (+1 to avoid 0 allocation when an empty surface comes in)
	surfacePoints = new double[3*numPoints+1];
	surfaceNormals = new double[3*numPoints+1];
	triangles = new int[3*numTriangles+1];

	///map in the data
	for(i = 0; i<numPoints; i++){
		surfacePoints[3*i+0] = ptsAndNorms[6*i+0];
		surfacePoints[3*i+1] = ptsAndNorms[6*i+1];
		surfacePoints[3*i+2] = ptsAndNorms[6*i+2];
		surfaceNormals[3*i+0] = ptsAndNorms[6*i+3];
		surfaceNormals[3*i+1] = ptsAndNorms[6*i+4];
		surfaceNormals[3*i+2] = ptsAndNorms[6*i+5];

		if(surfaceNormals[3*i+0] > 1.0 || surfaceNormals[3*i+0] < -1.0){
			printf("WTF: bad qhull norm\n");
			surfaceNormals[0] = 0;
		}
		if(surfaceNormals[3*i+1] > 1.0 || surfaceNormals[3*i+1] < -1.0){
			printf("WTF: bad qhull norm\n");
			surfaceNormals[1] = 0;
		}
		if(surfaceNormals[3*i+2] > 1.0 || surfaceNormals[3*i+2] < -1.0){
			printf("WTF: bad qhull norm\n");
			surfaceNormals[2] = 0;
		}
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
//		printf("TRIANGLE: %i %i %i\n", triangles[3*i+0], triangles[3*i+1], triangles[3*i+2] );
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
		//	printf("TEMPNORMAL: %f %f %f (point: %i)\n", tempNormal[0], tempNormal[1], tempNormal[2], triangles[3*i+0] );
		//	printf("CROSS: %f %f %f\n", cross[0], cross[1], cross[2]);
			
			CROSS(cross,v2,v1);
			//printf("CROSS: %f %f %f\n", cross[0], cross[1], cross[2]);
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
// Construction
//
// This is a secondary constructor, used for construction with existing data
// POints and norms are separate, as they would typically be in most programs.
SurfaceObject::SurfaceObject(int numPts, double * pts, double * norms, int numTris, int * tris, double * triangleNorms)
{
	int i = 0;
	numPoints = numPts;
	numTriangles = numTris;
	colors = NULL;
	highlights = NULL;
	edges = NULL;

	/////fill in the pts and averaged normals (+1 to avoid 0 allocation when an empty surface comes in)
	surfacePoints = new double[3*numPts+1];
	surfaceNormals = new double[3*numPts+1];
	for(i = 0; i<numPts; i++){
		surfacePoints[3*i+0] = pts[3*i+0];		surfacePoints[3*i+1] = pts[3*i+1];		surfacePoints[3*i+2] = pts[3*i+2];
		surfaceNormals[3*i+0] = norms[3*i+0];	surfaceNormals[3*i+1] = norms[3*i+1];	surfaceNormals[3*i+2] = norms[3*i+2];
	}

	///compute the centroid
	centroid = getCentroid();

	/////fill in the triangles (+1 to avoid 0 allocation when an empty surface comes in)
	triangles = new int[3*numTris+1];
	triangleNormals = new double[3*numTris+1];
	for(i = 0; i<numTris; i++){
		triangles[3*i+0] = tris[3*i+0];
		triangles[3*i+1] = tris[3*i+1];
		triangles[3*i+2] = tris[3*i+2];
		triangleNormals[3*i+0] = triangleNormals[3*i+0];
		triangleNormals[3*i+1] = triangleNormals[3*i+1];
		triangleNormals[3*i+2] = triangleNormals[3*i+2];
	}
}


///Assemble the surfaceObjects in the surfSet, return the result
SurfaceObject::SurfaceObject( set_t surfSet )
{
	int i = 0;
	int j = 0;
	int numSurfs = size_set(surfSet);

	//build an offset array to remember where all the points should be.
	int * pointOffset = new int[numSurfs];
	int * triangleOffset = new int[numSurfs];
	for(i = 0; i<numSurfs; i++){ pointOffset[i] = 0; triangleOffset[i] = 0; }
	for(i = 0; i<numSurfs-1; i++){
		SurfaceObject * surf = (SurfaceObject *) mapsto_set(surfSet, i);
		pointOffset[i+1] = pointOffset[i] + surf->numPoints;
		triangleOffset[i+1] = triangleOffset[i] + surf->numTriangles;
	}

	//sum the triangles
	numPoints = 0;
	numTriangles = 0;
	for(i = 0; i<numSurfs; i++){
		SurfaceObject * surf = (SurfaceObject *) mapsto_set(surfSet, i);
		numPoints += surf->numPoints;
		numTriangles += surf->numTriangles;
	}

	//now put everyuthing together.
	surfacePoints = new double[3*numPoints+1];
	surfaceNormals = new double[3*numPoints+1];
	triangles = new int[3*numTriangles+1];
	triangleNormals = new double[3*numTriangles+1];

	for(i = 0; i<numSurfs; i++){
		SurfaceObject * surf = (SurfaceObject *) mapsto_set(surfSet, i);
		for(j = 0; j<surf->numPoints; j++){
			surfacePoints[3*(pointOffset[i]+j)+0] = surf->surfacePoints[3*j+0];
			surfacePoints[3*(pointOffset[i]+j)+1] = surf->surfacePoints[3*j+1];
			surfacePoints[3*(pointOffset[i]+j)+2] = surf->surfacePoints[3*j+2];
			surfaceNormals[3*(pointOffset[i]+j)+0] = surf->surfaceNormals[3*j+0];
			surfaceNormals[3*(pointOffset[i]+j)+1] = surf->surfaceNormals[3*j+1];
			surfaceNormals[3*(pointOffset[i]+j)+2] = surf->surfaceNormals[3*j+2];
		}
		for(j = 0; j<surf->numTriangles; j++){
			triangles[3*(triangleOffset[i]+j)+0] = pointOffset[i] + surf->triangles[3*j+0];
			triangles[3*(triangleOffset[i]+j)+1] = pointOffset[i] + surf->triangles[3*j+1];
			triangles[3*(triangleOffset[i]+j)+2] = pointOffset[i] + surf->triangles[3*j+2];
			triangleNormals[3*(triangleOffset[i]+j)+0] = surf->triangleNormals[3*j+0];
			triangleNormals[3*(triangleOffset[i]+j)+1] = surf->triangleNormals[3*j+1];
			triangleNormals[3*(triangleOffset[i]+j)+2] = surf->triangleNormals[3*j+2];
		}
	}

	///compute the centroid
	centroid = getCentroid();
	colors = NULL;
	highlights = NULL;
	edges = NULL;

	delete[](pointOffset);
	delete[](triangleOffset);
}




//////////////////////////////////////////////////////////////////////
// Destruction
//
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
///returns the coords of the requested triangle
///returns null if out of range.
///notice that the triangle must be formatted in this way in the array (3 then 3 then 3)
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


///////////////////////////////////////////////////////
///Computes the centroid based on the point positions
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
///identifies edges which are not surrounded by triangles - locations where the surface
///is non manifold.
///This turns the edges variable non-null.
//
//NOTE: THIS FUNCTION IS CURRENTLY BROKEN
//
//
void SurfaceObject::identifyNonManifoldEdges(void)
{
	int i = 0;
	int j = 0;
	
	///create a set of sets to hold non-manifold edges
	set_t nonManifold = alloc_set(SP_MAP);
		
	///Create a set of sets to remember which triangles associate with which vertex.  Allocate that here.
	set_t verticesMappingToTris = alloc_set(SP_MAP);
	for(i = 0; i<numPoints; i++){
		set_t tris = alloc_set(0);
		verticesMappingToTris = associate_set(verticesMappingToTris, i, tris);
	}

	///associate each vertex with a list of triangles it is part of.
	///we do this by interating through all triangles
	set_t tris;
	for(i = 0; i<numTriangles; i++){
		tris = (set_t) mapsto_set(verticesMappingToTris, triangles[3*i+0]);
		tris = put_set(tris, i);
		verticesMappingToTris = associate_set(verticesMappingToTris, triangles[3*i+0], tris);

		tris = (set_t) mapsto_set(verticesMappingToTris, triangles[3*i+1]);
		tris = put_set(tris, i);
		verticesMappingToTris = associate_set(verticesMappingToTris, triangles[3*i+1], tris);
		
		tris = (set_t) mapsto_set(verticesMappingToTris, triangles[3*i+2]);
		tris = put_set(tris, i);
		verticesMappingToTris = associate_set(verticesMappingToTris, triangles[3*i+2], tris);
	}

	int min = 100000;
	int max = 0;
	int counter = 0;
	for(i = 0; i<size_set(verticesMappingToTris); i++){
		tris = (set_t) mapsto_set(verticesMappingToTris, verticesMappingToTris[i]);
		if(min > size_set(tris)){ 
			min = size_set(tris); 
			if(size_set(tris) == 0){
				counter++;
			}
		}
		if(max < size_set(tris)){ max = size_set(tris); }
	}

	printf("MAXIMUM NUMBER OF TRIANLGES ASSOCIATED WITH A point: %i\n", max);
	printf("MINIMUM NUMBER OF TRIANLGES ASSOCIATED WITH A point: %i  (total 0: %i)\n", min, counter);


	set_t set1;
	set_t set2;
	set_t tempSet;
	int numCommon;
	int * tri = new int[3];
	///iterate through all triangles, 
	for(i = 0; i<numTriangles; i++){
		tri[0] = triangles[3*i+0];
		tri[1] = triangles[3*i+1];
		tri[2] = triangles[3*i+2];

		///for each edge on the triangle, verify that exactly two 
		///distinct triangles are part of that edge.  otherwise it is
		///a non-manifold edge.
		for(j = 0; j<3; j++){
			//set the two vertices of the edge we are working on
			int t1 = tri[(j+0)%3];
			int t2 = tri[(j+1)%3];
			//identify the set of triangles adjacent to each endpoint
			set1 = (set_t) mapsto_set(verticesMappingToTris, t1);
			set2 = (set_t) mapsto_set(verticesMappingToTris, t2);
			///find out how many triangles are common between these endpoints
			numCommon = countNumCommon(set1, set2);

			//if numCommon == 1, then this is an edge of a nonmanifold surface.  Store it.
			if(numCommon == 1){
				if(!contains_set(nonManifold, t1)){ tempSet = alloc_set(0); }
				else{ tempSet = (set_t) mapsto_set(nonManifold, t1); }
				tempSet = put_set(tempSet, t2);
				nonManifold = associate_set(nonManifold, t1, tempSet);

				if(!contains_set(nonManifold, t2)){ tempSet = alloc_set(0); }
				else{ tempSet = (set_t) mapsto_set(nonManifold, t2); }
				tempSet = put_set(tempSet, t1);
				nonManifold = associate_set(nonManifold, t2, tempSet);
//				printf("Edge: %i %i\n", t1, t2);
			}

			//if numCommon is zero, then there is a serious issue with triangle topology.  Report it.
			if(numCommon == 0){
				printf("ERROR: Two vertices of the same triangle (vertex %i and %i in triangle %i) \n", t1, t2, i);
				printf("ERROR: are not members of at least one triangle, according to our data!!\n");
			}

			//if numCommon is greater than two, there is impossible surface occuring.  Report it.
			if(numCommon > 2){
				printf("ERROR: Vertices %i and %i appear to be adjacent to an edge which is part of \n", t1, t2);
				printf("ERROR: more than two triangles, creating a non-manifold surface! error!\n");
			}
			
			///if numCommon is equal to two, that is fine.  Report nothing.
		}

		///The end of this loop creates a list of non-manifold edges at the border of the surface, if it exists.
	}
	
	edges = nonManifold;
	
	delete[](tri);
	//free the data structure
	for(i = 0; i<numPoints; i++){ tris = (set_t) mapsto_set(verticesMappingToTris, i); free_set(tris); }
	free_set(verticesMappingToTris);
	
}





//////////////////////////////////////////////////////////////////////
///helper function tells you how int elements are identical between two sets.
///This probably should be in a library... but its too special purpose
int SurfaceObject::countNumCommon(set_t set1, set_t set2)
{
	int i = 0;
	int result = 0;

	int s1 = size_set(set1);		
	int s2 = size_set(set2);
	set_t tempSet1 = set1;		
	set_t tempSet2 = set2;

	if(s2 > s1){ 
		tempSet1 = set2;
		tempSet2 = set1;
		s1 = size_set(set2);
		s2 = size_set(set1);
	}
	
	///now tempSet1 is always at least as large as tempset2.
	for(i = 0; i<size_set(tempSet2); i++){
		if( contains_set(tempSet1, tempSet2[i]) ){
			result++;
		}
	}
	
	return result;
}





//////////////////////////////////////////////////////////////////////
///This function flips the normals of the surfaceObject backwards 
///so that we can treat it as a negative volume
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








//////////////////////////////////////////////////////////////////////
///adds color data to the object
void SurfaceObject::addColors(double * c)
{
	if(c == NULL){return;}
	
	colors = new double[3*numPoints];
	int i = 0;
	for(i = 0; i<3*numPoints; i++){
		colors[i] = c[i];
	}
}






///////////////////////////////////////////////////////
///Adds the geometry of another SurfaceObject to this SurfaceObject.
///ASSUMES THAT THIS SURFACEOBJECT DOES NOT COLLIDE WITH THE OTHER SURFACE OBJECT.
///DOES NOT CHECK FOR COLLISION, SO MUST BE DONE BEFORE THIS FUNCTION IS CALLED.
void SurfaceObject::addObject(SurfaceObject * obj)
{
	///Instantiate new copies of all the data in this class.
	int newNumPoints = numPoints + obj->numPoints;
//	printf("numNewPts: %i\n", newNumPoints);
	double * newSurfacePoints = new double[3*newNumPoints];
	double * newSurfaceNormals = new double[3*newNumPoints];
	int newNumTriangles = numTriangles + obj->numTriangles;
	int * newTriangles = new int[3*newNumTriangles];
	double * newTriangleNormals = new double[3*newNumTriangles];
	
//	printf("numPoints: %i  numTriangles: %i\n", numPoints, numTriangles);
//	printf("obj->numPoints: %i  obj->numTriangles: %i\n", obj->numPoints, obj->numTriangles);
//	printf("newNumPoints: %i  newNumTriangles: %i\n", newNumPoints, newNumTriangles);

	///Now copy over the old data and the expanded data for the points/normals
	int i = 0;
	int counter = 0;
	for(i = 0; i<numPoints; i++){
		newSurfacePoints[3*counter+0] = surfacePoints[3*i+0];
		newSurfacePoints[3*counter+1] = surfacePoints[3*i+1];
		newSurfacePoints[3*counter+2] = surfacePoints[3*i+2];
		newSurfaceNormals[3*counter+0] = surfaceNormals[3*i+0];
		newSurfaceNormals[3*counter+1] = surfaceNormals[3*i+1];
		newSurfaceNormals[3*counter+2] = surfaceNormals[3*i+2];
		counter++;
	}
	int newZeroCount = counter;
	for(i = 0; i<obj->numPoints; i++){
		newSurfacePoints[3*counter+0] = obj->surfacePoints[3*i+0];
		newSurfacePoints[3*counter+1] = obj->surfacePoints[3*i+1];
		newSurfacePoints[3*counter+2] = obj->surfacePoints[3*i+2];
		newSurfaceNormals[3*counter+0] = obj->surfaceNormals[3*i+0];
		newSurfaceNormals[3*counter+1] = obj->surfaceNormals[3*i+1];
		newSurfaceNormals[3*counter+2] = obj->surfaceNormals[3*i+2];
		counter++;
	}

	///Now copy over the old data and the expanded data for the triangles/normals
	counter = 0;
	for(i = 0; i<numTriangles; i++){
		newTriangles[3*counter+0] = triangles[3*i+0];
		newTriangles[3*counter+1] = triangles[3*i+1];
		newTriangles[3*counter+2] = triangles[3*i+2];
		newTriangleNormals[3*counter+0] = triangleNormals[3*i+0];
		newTriangleNormals[3*counter+1] = triangleNormals[3*i+1];
		newTriangleNormals[3*counter+2] = triangleNormals[3*i+2];
		counter++;
	}
	for(i = 0; i<obj->numTriangles; i++){
		newTriangles[3*counter+0] = newZeroCount + obj->triangles[3*i+0];
		newTriangles[3*counter+1] = newZeroCount + obj->triangles[3*i+1];
		newTriangles[3*counter+2] = newZeroCount + obj->triangles[3*i+2];
		newTriangleNormals[3*counter+0] = obj->triangleNormals[3*i+0];
		newTriangleNormals[3*counter+1] = obj->triangleNormals[3*i+1];
		newTriangleNormals[3*counter+2] = obj->triangleNormals[3*i+2];
		counter++;
	}

	///we delete these if they arent already NULL.
	if(highlights != NULL){
		free_set(highlights);  highlights = NULL;
	}
	if(edges != NULL){
		free_set(edges);  edges = NULL;
	}
	if(colors != NULL){
		delete[](colors);  colors = NULL;
	}

	///replace the existing data.	
	delete[](surfacePoints);
	surfacePoints = newSurfacePoints;

	delete[](surfaceNormals);
	surfaceNormals = newSurfaceNormals;

	delete[](triangles);
	triangles = newTriangles;

	delete[](triangleNormals);
	triangleNormals = newTriangleNormals;

	numPoints = newNumPoints;
	numTriangles = newNumTriangles;

	///compute the centroids.
	delete[](centroid);
	centroid = getCentroid();


}


///////////////////////////////////////////////////////
///Copies this surfaceObject into another surfaceObject
///If indices is non-NULL, returns a subset of the surface
///using only the points that are the indices of this surface.
SurfaceObject * SurfaceObject::copy(set_t indices)
{
	int i = 0;
	
	if(indices == NULL){
		SurfaceObject * result = new SurfaceObject(numPoints, surfacePoints, surfaceNormals, numTriangles, triangles, triangleNormals);
		return result;
	}

	//first gather the point data that is relevant.
	int newNumPoints = size_set(indices);
	double * newSurfacePoints = new double[3*newNumPoints];
	double * newSurfaceNormals = new double[3*newNumPoints];
	for(i = 0; i<newNumPoints; i++){
		int thisIndex = indices[i];
		newSurfacePoints[3*i+0] = surfacePoints[3*thisIndex+0];
		newSurfacePoints[3*i+1] = surfacePoints[3*thisIndex+1];
		newSurfacePoints[3*i+2] = surfacePoints[3*thisIndex+2];
		newSurfaceNormals[3*i+0] = surfaceNormals[3*thisIndex+0];
		newSurfaceNormals[3*i+1] = surfaceNormals[3*thisIndex+1];
		newSurfaceNormals[3*i+2] = surfaceNormals[3*thisIndex+2];
	}

	//now find the relevant triangles 
	set_t newTris = alloc_set(SP_MAP);
	set_t newNorms = alloc_set(SP_MAP);
	for(i = 0; i<numTriangles; i++){
		int t0 = triangles[3*i+0];
		int t1 = triangles[3*i+1];
		int t2 = triangles[3*i+2];

		bool t0Test = contains_set(indices, t0);
		bool t1Test = contains_set(indices, t1);
		bool t2Test = contains_set(indices, t2);

		///if the triangle is actually in the surface, then process it.
		if(t0Test && t1Test && t2Test){
			int * tris = new int[3];
			tris[0] = index_of_set(indices, t0);
			tris[1] = index_of_set(indices, t1);
			tris[2] = index_of_set(indices, t2);
			newTris = associate_set(newTris, size_set(newTris), tris);

			double * norm = new double[3];
			norm[0] = triangleNormals[3*i+0];
			norm[1] = triangleNormals[3*i+1];
			norm[2] = triangleNormals[3*i+2];
			newNorms = associate_set(newNorms, size_set(newNorms), norm);
		}
	}

	//map the relevant triangles into an acceptable array
	int newNumTriangles = size_set(newTris);
	int * newTriangles = new int[3*size_set(newTris)];
	double * newNormals = new double[3*size_set(newNorms)];
	for(i = 0; i<size_set(newTris); i++){
		int * tris = (int *) mapsto_set(newTris, i);
		newTriangles[3*i+0] = tris[0];
		newTriangles[3*i+1] = tris[1];
		newTriangles[3*i+2] = tris[2];

		double * norm = (double *) mapsto_set(newNorms, i);
		newNormals[3*i+0] = norm[0];
		newNormals[3*i+1] = norm[1];
		newNormals[3*i+2] = norm[2];
		
		delete[](tris);	//we dont need these anymore
		delete[](norm);
	}

	SurfaceObject * result = new SurfaceObject(newNumPoints, newSurfacePoints, newSurfaceNormals, newNumTriangles, newTriangles, newNormals);

	delete[](newTriangles);
	delete[](newNormals);
	free_set(newTris);
	free_set(newNorms);
	delete[](newSurfacePoints);
	delete[](newSurfaceNormals);

	return result;
}
///////////////////////////////////////////////////////




///////////////////////////////////////////////////////
///print out details
void SurfaceObject::toString()
{
	printf("Printing out Surface Object Data\n");

	int i = 0;
	for(i = 0; i<numPoints; i++){
		printf("Pt: (%f, %f, %f) Norm: (%f, %f, %f)\n", 
		surfacePoints[3*i+0], surfacePoints[3*i+1], surfacePoints[3*i+2], 
		surfaceNormals[3*i+0], surfaceNormals[3*i+1], surfaceNormals[3*i+2] );
	}

}


void SurfaceObject::printSummary()
{
	printf("SURFACE: [%i] points, [%i] triangles\n", numPoints, numTriangles);
		
}








