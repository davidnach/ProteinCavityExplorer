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
 * File: SurfaceOutput.cpp
 *       implementation for outputting triangle meshes
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

#include "SurfaceOutput.h"





////generateSURF:  This function generates a SURF file based on a surfaceObject and a set of
////                  point indices that reference a patch on the surface.  The output is a
////                  SURF file that contains only the points relevant to render the patch.
////                  If the set of points is NULL, the entire surface is returned.
////                  
////                  Inputs:
////                     a SurfaceObject *			describes the surface
////                     set_t points				describes which vertices we want (NULL for all)
////                     char * outputFileName         filename of the output surface file
////
////                  Outputs:
////                     Void - only a file is generated
////                  
////                  =========================================================
////                    OUTPUT FORMAT:
////                       -Anywhere in the file, empty lines and line starting with '#' are
////                       considered comments, and are ignored.
////                       -File has 2 sections: geometry and topology
////                  
////                       The geometry section is always first, and starts with:
////                       GEOMETRY: <int>
////                       where <int> is the number of points to be specified.
////                       then there are <int> lines as follows:
////                       <float> <float> <float> <float> <float> <float>
////                       which stand for x,y,z, and xnormal, ynormal, znormal.
////                  
////                       The topology section is always second, and starts with:
////                       TOPOLOGY: <int>
////                       where <int> specifies the number of triangles on the surface
////                       Then there are <int> lines as follows:
////                       <int> <int> <int>
////                       Each int stands for the 3 corners of the triangle
////                       and is an index into the array of points provided in the
////                       Geometry section.  I.e. the geometry section is indexed
////                       starting at zero, and ending at (size-1).
////                  =========================================================
////                  
void generateSURF(SurfaceObject * surf, set_t points, char * outputFileName)
{
	int i = 0;
	
	FILE * currentOutputFile = fopen(outputFileName, "w");
	
	int numVerts = surf->numPoints;
	if(points != NULL){
		numVerts = size_set(points);
	}

	if(points != NULL){
		printf("Generating SURF file [%s]: numVerts overall: %i, numVerts selected for output: %i\n", outputFileName, surf->numPoints, size_set(points));
	}
	else{
		printf("Generating SURF file [%s]: numVerts overall: %i, numVerts selected for output: %i\n", outputFileName, surf->numPoints, surf->numPoints);
	}

	fprintf(currentOutputFile, "#######################################################################\n" );
	fprintf(currentOutputFile, "## SURF file format 1.0                                              ##\n" );
	fprintf(currentOutputFile, "##   generated with VASP by Brian Chen, version %2.1f, 2014            ##\n", VERS_NUM );
	fprintf(currentOutputFile, "#######################################################################\n" );
	fprintf(currentOutputFile, "\n" );

	fprintf(currentOutputFile, "########################################\n" );
	fprintf(currentOutputFile, "###  ##   ##  ## ### #   #   #   ## # ##\n" );
	fprintf(currentOutputFile, "## #### ### ## #  #  # #### ## ## # # ##\n" );
	fprintf(currentOutputFile, "## #  #  ## ## # # # #  ### ##   ### ###\n" );
	fprintf(currentOutputFile, "## ## # ### ## # ### # #### ## # ### ###\n" );
	fprintf(currentOutputFile, "###  ##   ##  ## ### #   ## ## ## ## ###\n" );
	fprintf(currentOutputFile, "########################################\n" );
	fprintf(currentOutputFile, "GEOMETRY: %i\n", numVerts);


	///generate the Geometry
	int * map = NULL;
	int * rmap = NULL;
	///if points is null, jsut output the geometry
	if(points==NULL){
		for(i = 0; i<numVerts; i++){
			fprintf(currentOutputFile, "%f %f %f %f %f %f\n", 
				surf->surfacePoints[3*i+0], surf->surfacePoints[3*i+1], surf->surfacePoints[3*i+2], 
				surf->surfaceNormals[3*i+0], surf->surfaceNormals[3*i+1], surf->surfaceNormals[3*i+2] );
		}
	}
	//if points is not null, prepare forward and backward maps so that we cna quickly process the triangles.
	//then generate only the points used by this pocket.
	else{
		map = new int[size_set(points)];
		rmap = new int[surf->numPoints];
		for(i = 0; i<surf->numPoints; i++){ rmap[i] = -1; }
		for(i = 0; i<size_set(points); i++){
			map[i] = points[i];
			rmap[points[i]] = i;
		}
		for(i = 0; i<numVerts; i++){
			int tempVal = points[i];
			fprintf(currentOutputFile, "%f %f %f %f %f %f\n", 
				surf->surfacePoints[3*tempVal+0], surf->surfacePoints[3*tempVal+1], surf->surfacePoints[3*tempVal+2], 
				surf->surfaceNormals[3*tempVal+0], surf->surfaceNormals[3*tempVal+1], surf->surfaceNormals[3*tempVal+2] );
		}
	}

	///generate the Topology Header
	fprintf(currentOutputFile, "########################################\n" );
	fprintf(currentOutputFile, "##   ##  ##   ###  ## ####  ###  ## # ##\n" );
	fprintf(currentOutputFile, "### ## ## # ## # ## # ### ## # #### # ##\n" );
	fprintf(currentOutputFile, "### ## ## #   ## ## # ### ## # #  ## ###\n" );
	fprintf(currentOutputFile, "### ## ## # #### ## # ### ## # ## ## ###\n" );
	fprintf(currentOutputFile, "### ###  ## #####  ##   ##  ###  ### ###\n" );
	fprintf(currentOutputFile, "########################################\n" );


	///count the number of triangles:
	int numTriangles = 0;
	///If points is null, generate all triangles
	if(points==NULL){
		numTriangles = surf->numTriangles;
	}
	///If points is not null, count the triangles list for the triangles which are
	///entirely contained within points, and only output those.
	else{
		for(i = 0; i<surf->numTriangles; i++){
			bool test1 = contains_set(points, surf->triangles[3*i+0]);
			bool test2 = contains_set(points, surf->triangles[3*i+1]);
			bool test3 = contains_set(points, surf->triangles[3*i+2]);
			if(test1 && test2 && test3){
				numTriangles++;
			}
		}		
	}	
	///output the number of triangles.	
	fprintf(currentOutputFile, "TOPOLOGY: %i\n", numTriangles);

	//now actually output the triangles that we need.Use the reverse map.
	int point1Idx, point2Idx, point3Idx;
	///if points is null, just output everything.
	if(points==NULL){
		for(i = 0; i<numTriangles; i++){
			///get the triangle indices from surfTriangles
			point1Idx = surf->triangles[3*i+0];
			point2Idx = surf->triangles[3*i+1];
			point3Idx = surf->triangles[3*i+2];
			fprintf(currentOutputFile, "%i %i %i\n", point1Idx, point2Idx, point3Idx);
		}
	}
	//If poitns is non-null, output only the ones we care about.  Use the reverse map.
	else{
		for(i = 0; i<surf->numTriangles; i++){
			bool test1 = contains_set(points, surf->triangles[3*i+0]);
			bool test2 = contains_set(points, surf->triangles[3*i+1]);
			bool test3 = contains_set(points, surf->triangles[3*i+2]);
			if(test1 && test2 && test3){
				fprintf(currentOutputFile, "%i %i %i\n", rmap[surf->triangles[3*i+0]], rmap[surf->triangles[3*i+1]], rmap[surf->triangles[3*i+2]]);
			}
		}
	}

	if(surf->colors != NULL){	
		fprintf(currentOutputFile, "\n" );
		fprintf(currentOutputFile, "########################################\n" );
		fprintf(currentOutputFile, "###  ##  ## ####  ##  ###   ############\n" );
		fprintf(currentOutputFile, "## ### ## # ### ## # # # ###############\n" );
		fprintf(currentOutputFile, "## ### ## # ### ## #  ###  #############\n" );
		fprintf(currentOutputFile, "## ### ## # ### ## # # #### ############\n" );
		fprintf(currentOutputFile, "###  ##  ##   ##  ## # #   #############\n" );
		fprintf(currentOutputFile, "########################################\n" );
		fprintf(currentOutputFile, "COLORS: \n" );

		for(i = 0; i<surf->numPoints; i++){
//			sprintf(tempString, "%f %f %f\n", surf->colors[3*i+0], surf->colors[3*i+1], surf->colors[3*i+2] );
//			( *currentOutputFile ) << tempString;
			fprintf(currentOutputFile, "%f %f %f\n", surf->colors[3*i+0], surf->colors[3*i+1], surf->colors[3*i+2] );
		}
	}

	///now output the footer.
	fprintf(currentOutputFile, "\n" );
	fprintf(currentOutputFile, "#######################################################################\n" );
	fprintf(currentOutputFile, "## SURF file format 1.0                                              ##\n" );
	fprintf(currentOutputFile, "##   generated with VASP by Brian Chen, version %2.1f, 2014            ##\n", VERS_NUM );
	fprintf(currentOutputFile, "#######################################################################\n" );

	///done outputting.  Close file.
	fflush(currentOutputFile);
	fclose(currentOutputFile);
	
	printf("Output SURF File [%s] complete.\n", outputFileName);

	if(map != NULL){
		delete[](map);
	}
	if(rmap != NULL){
		delete[](rmap);
	}
}



























//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////
////Helper function Interpolates between two points with two charges, tells where the threshold point is.
double * pInterpolate(double c1, double * p1, double c2, double * p2, double zeroValue)
{
//	if( (c1 != c1) || (c2 != c2) ){
//		printf("WTF! nan came back!\n");
//		exit(1);
//	}

	double * thisPoint = new double[3];
	double d1 = c1;
	double d2 = c2;

	if( d1 != d1 ){ d1 = 0.0; }	
	if( d2 != d2 ){ d2 = 0.0; }

//	printf("%f %f %f (%f)   %f %f %f (%f)    %f\n", p1[0], p1[1], p1[2], c1, p2[0], p2[1], p2[2], c2, zeroValue );
//	exit(1);

	double diff0 = p1[0] - p2[0];
	double diff1 = p1[1] - p2[1];
	double diff2 = p1[2] - p2[2];
	
	if(diff0 == 0){ thisPoint[0] = p1[0]; } 
	else{ 
		if(diff0 < 0){ 
			thisPoint[0] = p1[0] + ((zeroValue-d1)/(d1-d2))*(p1[0]-p2[0]); 
		}
		else{ 
			thisPoint[0] = p2[0] + ((zeroValue-d2)/(d2-d1))*(p2[0]-p1[0]); 
		}
	}
	if(diff1 == 0){ thisPoint[1] = p1[1]; } 
	else{ 
		if(diff1 < 0){ 
			thisPoint[1] = p1[1] + ((zeroValue-d1)/(d1-d2))*(p1[1]-p2[1]); 
		}
		else{ 
			thisPoint[1] = p2[1] + ((zeroValue-d2)/(d2-d1))*(p2[1]-p1[1]);
		}
	}
	if(diff2 == 0){ thisPoint[2] = p1[2]; } 
	else{ 
		if(diff2 < 0){
			thisPoint[2] = p1[2] + ((zeroValue-d1)/(d1-d2))*(p1[2]-p2[2]);
		}
		else{
			thisPoint[2] = p2[2] + ((zeroValue-d2)/(d2-d1))*(p2[2]-p1[2]);
		}
	}

	return thisPoint;	
}

bool pCheck(double eVal, double threshold, int rangeMode)
{
	bool result = false;
	if(rangeMode == -1 && eVal < threshold){ result = true; }
	if(rangeMode ==  1 && eVal > threshold){ result = true; }

	return result;
}


//////////////////////////////////////////////////////////
////Surface Generation based on values in INSIGHT grid
///
///  if rangeMode == -1, then generate the volume surrounding the region below the threshold,
///  if rangeMode == +1, then generate the volume surrounding the region above the threshold,
SurfaceObject * generatePhimapSurface(char * pdbFile, char * insightGrid, double thresh, int rangeMode, double res)
{
	
	//error checking
	if(rangeMode !=-1 && rangeMode !=1){ printf("ERROR: rangeMode value [%i] makes no sense.  Exitting.\n", rangeMode); exit(1); }

	//test if the pdb file exists.
	FILE * testFile = fopen(pdbFile, "r");
	if(testFile == NULL){ printf("ERROR: pdb File [%s] does not exist.  Exitting.\n", pdbFile); exit(1); }
	else{ fclose(testFile); }

	//test if the insight file exists.
	testFile = fopen(insightGrid, "r");
	if(testFile == NULL){ printf("ERROR: INSIGHT file [%s] does not exist.  Exitting.\n", insightGrid); exit(1); }
	else{ fclose(testFile); }
	
	///generate requisite data structures
//	double hashRes = 2.0;
//	if(res < 2.0){
//		hashRes = (res+2.0)/2;
//	}

//	TriangleLatticeHash * hash = new TriangleLatticeHash(surf, hashRes, "SurfaceToSimplify");
	CubeTable * table = new CubeTable();
//	table->readFile();

	///initialize counters
	int i = 0;	int j = 0;	int k = 0;	int l = 0;

	///get boundary values
	double xneg = HUGE_VAL;		double xpos = -HUGE_VAL;
	double yneg = HUGE_VAL;		double ypos = -HUGE_VAL;
	double zneg = HUGE_VAL;		double zpos = -HUGE_VAL;

//	////////////////////////////////////////////////////
//	////parse the potentials data.  It is stored in pmap
//	printf("Parsing potentials from INSIGHT grid file...\n"); fflush(stdout);
//	PDBFile mypdbFile(pdbFile, allatoms);
//	Structure mystruct;
//	mypdbFile >> mystruct;
//	Phimap * pmap = new Phimap(&mystruct, allatoms);
//	pmap->LoadInsight(insightGrid);
//	////parse the potentials data.  It is stored in pmap
//	////////////////////////////////////////////////////

	////////////////////////////////////////////////////
	////parse the potentials data.  It is stored in pmap
	printf("Parsing potentials from INSIGHT grid file...\n"); fflush(stdout);
	phimapGrid * pmap = new phimapGrid(insightGrid);
	////parse the potentials data.  It is stored in pmap
	////////////////////////////////////////////////////

	xneg = pmap->origin[0];
	yneg = pmap->origin[1];
	zneg = pmap->origin[2];
	xpos = pmap->extent[0];
	ypos = pmap->extent[1];
	zpos = pmap->extent[2];
	
	////////////////////////////////////////
	printf("PHIMAP DIMENSIONS:  xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	////////////////////////////////////////

	///negative padding to force a fit.
	xneg += 5;
	yneg += 5;
	zneg += 5;
	xpos -= 5;
	ypos -= 5;
	zpos -= 5;

	//padding on low side.
	xneg -= res;	yneg -= res;	zneg -= res;
	///integer iteration
	xneg = floor(xneg);		yneg = floor(yneg);		zneg = floor(zneg);
	xpos = ceil(xpos);		ypos = ceil(ypos);		zpos = ceil(zpos);

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
	printf("Computing Simplification...\n"); fflush(stdout);
	int pointCounter = 0;
	set_t points = alloc_set(SP_MAP);
	for(i = 0; i<xdim; i++){
		for(j = 0; j<ydim; j++){
			for(k = 0; k<zdim; k++){
				stdoutProgressBar(progress_cube, progress_max);

				///Allocate Stuff we will use.
				//int oldPtCounter = pointCounter;
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

				///store electrostatic potentials here.
				double * e = new double[8];
				e[0] = 0.0;	e[1] = 0.0;	e[2] = 0.0;	e[3] = 0.0;	e[4] = 0.0;	e[5] = 0.0;	e[6] = 0.0;	e[7] = 0.0;
				
				try{ e[0] = pmap->getPotential(c0[0], c0[1], c0[2]); if (e[0]!=e[0]){e[0]=0.0;} /* if(e[0]!=0){printf("<<%f>>\n", e[0]);}*/ } catch(...) { e[0] = 0.0; }
				try{ e[1] = pmap->getPotential(c1[0], c1[1], c1[2]); if (e[1]!=e[1]){e[1]=0.0;} /* if(e[1]!=0){printf("<<%f>>\n", e[1]);}*/ } catch(...) { e[1] = 0.0; }
				try{ e[2] = pmap->getPotential(c2[0], c2[1], c2[2]); if (e[2]!=e[2]){e[2]=0.0;} /* if(e[2]!=0){printf("<<%f>>\n", e[2]);}*/ } catch(...) { e[2] = 0.0; }
				try{ e[3] = pmap->getPotential(c3[0], c3[1], c3[2]); if (e[3]!=e[3]){e[3]=0.0;} /* if(e[3]!=0){printf("<<%f>>\n", e[3]);}*/ } catch(...) { e[3] = 0.0; }
				try{ e[4] = pmap->getPotential(c4[0], c4[1], c4[2]); if (e[4]!=e[4]){e[4]=0.0;} /* if(e[4]!=0){printf("<<%f>>\n", e[4]);}*/ } catch(...) { e[4] = 0.0; }
				try{ e[5] = pmap->getPotential(c5[0], c5[1], c5[2]); if (e[5]!=e[5]){e[5]=0.0;} /* if(e[5]!=0){printf("<<%f>>\n", e[5]);}*/ } catch(...) { e[5] = 0.0; }
				try{ e[6] = pmap->getPotential(c6[0], c6[1], c6[2]); if (e[6]!=e[6]){e[6]=0.0;} /* if(e[6]!=0){printf("<<%f>>\n", e[6]);}*/ } catch(...) { e[6] = 0.0; }
				try{ e[7] = pmap->getPotential(c7[0], c7[1], c7[2]); if (e[7]!=e[7]){e[7]=0.0;} /* if(e[7]!=0){printf("<<%f>>\n", e[7]);}*/ } catch(...) { e[7] = 0.0; }

				//////////////////////////////////////////////////////////////////////////////////////////////
				///SPECIAL POTENTIAL CONTOUR BORDER PROCESSING
				//////////////////////////////////////////////////////////////////////////////////////////////
				///Potential contours can be infinite.  Check to see if this is a border cube, and if it is,
				///Ensure that the outside border of the cube is outside the region.
				///
				///   3---7            ^
				///   |\  |\       \   |
				///   | 1---5       +y |
				///   | | | |        \ z+
				///   2-|-6 |         \|  
				///    \|  \|          0--x+-->    
				///     0---4
				///
				double borderVal;
				if(rangeMode==-1){ borderVal = thresh+.001; }
				if(rangeMode== 1){ borderVal = thresh-.001; }
				if( i == 0 ){ e[0]=borderVal; e[1]=borderVal; e[2]=borderVal; e[3]=borderVal; }  ///lowest Xdim cube
				if( j == 0 ){ e[0]=borderVal; e[1]=borderVal; e[4]=borderVal; e[5]=borderVal; }  ///lowest Ydim cube
				if( k == 0 ){ e[0]=borderVal; e[2]=borderVal; e[4]=borderVal; e[6]=borderVal; }  ///lowest Zdim cube
				if( i == xdim-1 ){ e[4]=borderVal; e[5]=borderVal; e[6]=borderVal; e[7]=borderVal; }  ///highest Xdim cube
				if( j == ydim-1 ){ e[2]=borderVal; e[3]=borderVal; e[6]=borderVal; e[7]=borderVal; }  ///highest Ydim cube
				if( k == zdim-1 ){ e[1]=borderVal; e[3]=borderVal; e[5]=borderVal; e[7]=borderVal; }  ///highest Zdim cube
				//////////////////////////////////////////////////////////////////////////////////////////////
				///SPECIAL POTENTIAL CONTOUR BORDER PROCESSING
				//////////////////////////////////////////////////////////////////////////////////////////////

				///Test if the points fall inside or outside the region.
				bool p0 = false, p1 = false, p2 = false, p3 = false, p4 = false, p5 = false, p6 = false, p7 = false;
				if( pCheck(e[0], thresh, rangeMode) ){ p0 = true; thesePts[0] = true;}  // printf("pt 0 is inside\n"); }
				if( pCheck(e[1], thresh, rangeMode) ){ p1 = true; thesePts[1] = true;}  // printf("pt 1 is inside\n"); }
				if( pCheck(e[2], thresh, rangeMode) ){ p2 = true; thesePts[2] = true;}  // printf("pt 2 is inside\n"); }
				if( pCheck(e[3], thresh, rangeMode) ){ p3 = true; thesePts[3] = true;}  // printf("pt 3 is inside\n"); }
				if( pCheck(e[4], thresh, rangeMode) ){ p4 = true; thesePts[4] = true;}  // printf("pt 4 is inside\n"); }
				if( pCheck(e[5], thresh, rangeMode) ){ p5 = true; thesePts[5] = true;}  // printf("pt 5 is inside\n"); }
				if( pCheck(e[6], thresh, rangeMode) ){ p6 = true; thesePts[6] = true;}  // printf("pt 6 is inside\n"); }
				if( pCheck(e[7], thresh, rangeMode) ){ p7 = true; thesePts[7] = true;}  // printf("pt 7 is inside\n"); }

				///SET THE EDGES BASED ON THE INSIDE/OUTSIDE PARITY
				///EDGE 0//////////////////////////////////////////////////////////////////////////////////////
				if((p0 && !p1) || (!p0 && p1)){
					thisPoint = pInterpolate(e[0], c0, e[1], c1, thresh);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[0] = pointCounter; pointCounter++;
				}
				///EDGE 1//////////////////////////////////////////////////////////////////////////////////////
				if((p0 && !p2) || (!p0 && p2)){
					thisPoint = pInterpolate(e[0], c0, e[2], c2, thresh);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[1] = pointCounter; pointCounter++;
				}
				///EDGE 4//////////////////////////////////////////////////////////////////////////////////////
				if((p0 && !p4) || (!p0 && p4)){
					thisPoint = pInterpolate(e[0], c0, e[4], c4, thresh);
					points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
					theseEdges[4] = pointCounter; pointCounter++;
				}
				////If you are on one of the high X side, get the high X edge
				if(i == xdim-1){
				///EDGE 8//////////////////////////////////////////////////////////////////////////////////////
					if((p4 && !p5) || (!p4 && p5)){
						thisPoint = pInterpolate(e[4], c4, e[5], c5, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[8] = pointCounter; pointCounter++;
					}
				///EDGE 9//////////////////////////////////////////////////////////////////////////////////////
					if((p4 && !p6) || (!p4 && p6)){
						thisPoint = pInterpolate(e[4], c4, e[6], c6, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[9] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Y side, get the high Y edge
				if(j == ydim-1){
				///EDGE 3//////////////////////////////////////////////////////////////////////////////////////
					if((p2 && !p3) || (!p2 && p3)){
						thisPoint = pInterpolate(e[2], c2, e[3], c3, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[3] = pointCounter; pointCounter++;
					}
				///EDGE 6//////////////////////////////////////////////////////////////////////////////////////
					if((p2 && !p6) || (!p2 && p6)){
						thisPoint = pInterpolate(e[2], c2, e[6], c6, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[6] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Z side, get the high Z edge
				if(k == zdim-1){
				///EDGE 2//////////////////////////////////////////////////////////////////////////////////////
					if((p1 && !p3) || (!p1 && p3)){
						thisPoint = pInterpolate(e[1], c1, e[3], c3, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[2] = pointCounter; pointCounter++;
					}
				///EDGE 5//////////////////////////////////////////////////////////////////////////////////////
					if((p1 && !p5) || (!p1 && p5)){
						thisPoint = pInterpolate(e[1], c1, e[5], c5, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[5] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high X and Y sides, get the high X and Y edge
				if(i == xdim-1 && j == ydim-1){
				///EDGE 11//////////////////////////////////////////////////////////////////////////////////////
					if((p6 && !p7) || (!p6 && p7)){
						thisPoint = pInterpolate(e[6], c6, e[7], c7, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[11] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high Y and Z sides, get the high Y and Z edge
				if(j == ydim-1 && k == zdim-1){
				///EDGE 7//////////////////////////////////////////////////////////////////////////////////////
					if((p3 && !p7) || (!p3 && p7)){
						thisPoint = pInterpolate(e[3], c3, e[7], c7, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[7] = pointCounter; pointCounter++;
					}
				}
				////If you are on one of the high X and Z sides, get the high X and Z edge
				if(k == zdim-1 && i == xdim-1){
				///EDGE 10//////////////////////////////////////////////////////////////////////////////////////
					if((p5 && !p7) || (!p5 && p7)){
						thisPoint = pInterpolate(e[5], c5, e[7], c7, thresh);
						points = associate_set(points, pointCounter, (ptr_t) thisPoint); 
						theseEdges[10] = pointCounter; pointCounter++;
					}
				}
				
				delete[](e);

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
	printf("Generating Triangles and Normals.\n"); fflush(stdout);
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
			//double * thisVector = (double *) mapsto_set(points, thisPt);  ///this gets the coordinates of the point
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

	////STEP 3a: Generate Electrostatic Coloring
	double * colors = new double[3*numPts];
	double * colorVal = new double[3];
	
	//double thresh, int rangeMode, double res)
	
	if(rangeMode == -1){ colorVal[0] = 1; colorVal[1] = 0; colorVal[2] = 0;  }
	if(rangeMode == +1){ colorVal[0] = 0; colorVal[1] = 0; colorVal[2] = 1;  }
		
	for(i = 0; i<numPts; i++){
		colors[3*i+0] = colorVal[0];
		colors[3*i+1] = colorVal[1];
		colors[3*i+2] = colorVal[2];
	}
	delete[](colorVal);


	////STEP 4: Clean up all the memory allocations
	///////////////////////////////////////////////////////////////////////////////
	for(i = 0; i<size_set(points); i++){ double * vector = (double *) mapsto_set(points, i); delete[](vector); }
	for(i = 0; i<size_set(cubes); i++){ int * cube = (int *) mapsto_set(cubes, i); delete[](cube); }
	for(i = 0; i<triListSize; i++){ delete[](triList[i]); delete[](normList[i]); }
	free_set(points);
	delete[](triList);
	delete[](normList);
	delete(triSet);
	///////////////////////////////////////////////////////////////////////////////
	////STEP 4: Clean up all the memory allocations
	
	////STEP 5: Generate the output:	
	///////////////////////////////////////////////////////////////////////////////
	SurfaceObject * results = new SurfaceObject (numPts, pts, norms, numTris, tris, triangleNorms);
	results->addColors(colors);

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
	delete[](colors);
	///////////////////////////////////////////////////////////////////////////////
	////Step 6: Clear the output data:


	/////Return the result
	///////////////////////////////////////////////////////////////////////////////
	return results;





}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////











