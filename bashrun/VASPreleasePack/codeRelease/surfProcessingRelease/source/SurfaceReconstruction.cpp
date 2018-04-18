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
 * File: SurfaceReconstruction.h
 *       measures volume inside the surface
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


#include "SurfaceReconstruction.h"





//////////////////////
////Surveyor's Formula
////-- This volume measurement algorithm introduced to satisfy reviewers for PLOS Comp BIol.
void surveyorsFormulaExecution(char * objFileName)
{
	printf("Parsing File [%s]...\n", objFileName);
	SurfaceObject * obj = parseGeometryFile(objFileName); 
	printf("Running Surveyor's Formula Volume Computation...\n");
	double vol = computeVolumeSurveyorsFormula(obj);
	printf("Final Volume [%s]: %f\n", objFileName, vol);
}

//////////////////////
////Surveyor's Formula
////-- This volume measurement algorithm introduced to satisfy reviewers for PLOS Comp BIol.
double computeVolumeSurveyorsFormula(SurfaceObject * obj)
{
	int i = 0;
	double volume = 0;
	double * v0 = new double[3];	///these three vertices are used to store the triangle vertices.
	double * v1 = new double[3];	/// so that we dont have to do lots of news and deletes in the loop.
	double * v2 = new double[3];
	double * triNorm = new double[3];	///this vector stores the triangle normal as we use it.
	double * tCent = new double[3];
	double * nVec = new double[3]; //normalizedVector;

	///first define the centriod of the volume.
	double * c = obj->getCentroid(); 
	c[0] = obj->centroid[0];
	c[1] = obj->centroid[1];
	c[2] = obj->centroid[2];

//	c[0] += 1000;
//	c[1] += 1000;
//	c[2] += 1000;
	
	///second, iterate through all the triangles of the object.
	for(i = 0; i<obj->numTriangles; i++){
		///record the indices for this triangle.
		int t0 = obj->triangles[3*i+0];
		int t1 = obj->triangles[3*i+1];
		int t2 = obj->triangles[3*i+2];

		///range check.
		if( t0>obj->numPoints || t1>obj->numPoints || t2>obj->numPoints ){
			printf("ERROR: index %i, %i, or %i is out of range (%i)! exitting.\n", t0, t1, t2, obj->numPoints);
		}

		//record the surface points		
		v0[0] = obj->surfacePoints[3*t0+0];	v0[1] = obj->surfacePoints[3*t0+1];	v0[2] = obj->surfacePoints[3*t0+2];
		v1[0] = obj->surfacePoints[3*t1+0];	v1[1] = obj->surfacePoints[3*t1+1];	v1[2] = obj->surfacePoints[3*t1+2];
		v2[0] = obj->surfacePoints[3*t2+0];	v2[1] = obj->surfacePoints[3*t2+1];	v2[2] = obj->surfacePoints[3*t2+2];
		
		//record the triangle norm
		triNorm[0] = obj->triangleNormals[3*i+0];
		triNorm[1] = obj->triangleNormals[3*i+1];
		triNorm[2] = obj->triangleNormals[3*i+2];
		
		//compute the vector from the centroid of the triangle to the object centroid.
		tCent[0] = (v0[0] + v1[0] + v2[0])/3 - c[0];
		tCent[1] = (v0[1] + v1[1] + v2[1])/3 - c[1];
		tCent[2] = (v0[2] + v1[2] + v2[2])/3 - c[2];

		///normalize the vector
		normalizeVector(tCent, nVec);

		///now compute the dot product with the current vector.
		///If dProd is positive, then the two vectors face in the 
		///same dir, negative opposite dir, 0, perpendicular
		double dProd = dotProd(triNorm, nVec);
		if( dProd>0 ){ dProd = 1; }
		else if( dProd<0 ){ dProd = -1; }
		else{ dProd = 0; }

		///now compute the volume of the tetrahedron filled by the triangle and the centroid.
		double thisVol = computeTetrahedralVolume( v0, v1, v2, c);

		///add the volume (or subtract the volume) from the total.
		if(dProd != 0){
			//printf("Volume: %f  Volume Change: %f\n", volume, dProd*thisVol);
			volume += dProd*thisVol;
		}
	}
	
	
	delete[](c);
	delete[](v0);
	delete[](v1);
	delete[](v2);
	delete[](triNorm);
	delete[](tCent);
	delete[](nVec);
	
	return volume;
}
////Surveyor's Formula
////-- This volume measurement algorithm introduced to satisfy reviewers for PLOS Comp BIol.
//////////////////////








//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///input is four points in 3d space, implemented volume on the method by Niccolò Fontana Tartaglia, which
///is a generalization of Heron's area method:
//                 [  0     1    1     1     1    ]
//                 [  1     0   d12^2 d13^2 d14^2 ]
// 288 * v^2 = det [  1   d21^2  0    d23^2 d24^2 ]
//                 [  1   d31^2 d32^2  0    d34^2 ]
//                 [  1   d41^2 d42^2 d43^2  0    ]
//
//OPTIMIZATION: implement this so that it runs the 4x4 matrices at a more optimal line than the first.
//              (not coded yet)
//

double computeTetrahedralVolume( double * p1, double * p2, double * p3, double * p4 )
{
	/////////////////////////////////
	double d12 = vectorSize( p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2] );
	double d13 = vectorSize( p1[0]-p3[0], p1[1]-p3[1], p1[2]-p3[2] );
	double d14 = vectorSize( p1[0]-p4[0], p1[1]-p4[1], p1[2]-p4[2] );
	
	double d23 = vectorSize( p2[0]-p3[0], p2[1]-p3[1], p2[2]-p3[2] );
	double d24 = vectorSize( p2[0]-p4[0], p2[1]-p4[1], p2[2]-p4[2] );
	
	double d34 = vectorSize( p3[0]-p4[0], p3[1]-p4[1], p3[2]-p4[2] );
	/////////////////////////////////
	double s12 = d12*d12;
	double s13 = d13*d13;
	double s14 = d14*d14;

	double s23 = d23*d23;
	double s24 = d24*d24;

	double s34 = d34*d34;
	/////////////////////////////////
	double * m = new double[25];
	m[ 0] = 0;		m[ 5] = 1;		m[10] = 1;		m[15] = 1;		m[20] = 1;
	m[ 1] = 1;		m[ 6] = 0;		m[11] = s12;	m[16] = s13;	m[21] = s14;
	m[ 2] = 1;		m[ 7] = s12;	m[12] = 0;		m[17] = s23;	m[22] = s24;
	m[ 3] = 1;		m[ 8] = s13;	m[13] = s23;	m[18] = 0;		m[23] = s34;
	m[ 4] = 1;		m[ 9] = s14;	m[14] = s24;	m[19] = s34;	m[24] = 0;
	/////////////////////////////////
	double det = determinant5x5(m);

	delete[](m);

	double result = sqrt(det/288);
	
	if (result != result){
		printf("WARNING: ComputeTetrahedralVolume generated a NAN tetrahedron!\n");
		result = 0;
	}
	
	
	return result;
}
