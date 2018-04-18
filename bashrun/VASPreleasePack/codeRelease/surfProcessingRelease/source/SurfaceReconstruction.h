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
 *       interface for measuring volume inside the surface
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

#ifndef SURFACE_RECONSTRUCTION
#define SURFACE_RECONSTRUCTION

#include "StdAfx.h"

/*
	///Uses union find to identify interior cavities (disconnected surfaces)
	///and eliminate them.  Expects that the structure with the most triangles
	///is the exterior surface.  (this is a short cut)
	int * eliminateInteriorCavities(double * pts, int numPts, int * triangles, int numberTriangles, int * numTriangles);


	///this union find groups together all connected patches, and returns them in a set of sets
	set_t findDisparatePatches(double * inputPts, int inputNumPts, int * inputTriangles, int numberTriangles);


	///separate fragments stuff here
	void sepFragsExecution(char * surfFile, char * outputHeader);
	void surfaceTopologyHelper(set_t outsideSets, int * depths, int * flags, int currentSet, int currentDepth, int currentFlag);
	bool surfInsideSurf(SurfaceObject * s1, TriangleLatticeHash * h1, SurfaceObject * s2);

	///this function surveys the areas of triangles in the surface, and outputs the results in a histogram.
	///returns total area.  No output if numBars == 0.
	double surveyTriangleAreas(SurfaceObject * obj, int numBars);

	///this function surveys the edge Lengths of triangles in the surface, and outputs the results in a histogram. (all edges appear twice)
	void surveyTriangleEdges(SurfaceObject * obj, int numBars);
*/
	//////////////////////
	////Surveyor's Formula
	////-- This volume measurement algorithm introduced to satisfy reviewers for PLOS Comp BIol.
	double computeVolumeSurveyorsFormula(SurfaceObject * obj);
	void surveyorsFormulaExecution(char * objFileName);


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
	double computeTetrahedralVolume( double * p1, double * p2, double * p3, double * p4 );
	
	
	
	//////////////////////
	////Compute residue level Solvent Accessibility, output to stdout.
	////output works like this: FILENAME   RES#   Solv.Ass.Area   TotalSurfaceArea
	////IMPORTANT NOTE: WILL NOT WORK IF THERE ARE MULTIPLE CHAINS IN THIS PDB FILE
	////THIS ASSUMES THAT ALL RESIDUE NUMBERS ARE UNIQUE.  IF TEHY ARE NOT, THEN ALL
	////AMINO ACIDS WITH THE SAME ACCESSIBILITY WILL HAVE SURFACE AREA MEASURED AT ONCE
//	void computeSolventAccessibility( char * pFile );
	


#endif



