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
 * File: GeometryParser.cpp
 *       SURF file parser
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

#include "GeometryParser.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//
//NOTES: this file will need to be updated for more general parsing.
//In particular, the Geometry, Topology, Colors, etc sections must be allowed to have arbitrary orders.
//
//
//
//
SurfaceObject * parseGeometryFile(char * fileName)
{
	FILE * surfFile = fopen(fileName, "r");
	if(surfFile == NULL){
		printf("ERROR: Surface file [%s] not found.  Exitting.\n", fileName);
		exit(1);
	}
	else{
		printf("Parsing Surface Geometry File: [%s]\n", fileName);
	}

	int numPoints = getNumberOfPoints(surfFile);
	double * VectorsAndNormals = getSurfaceGeometryAndNormals(surfFile, numPoints);
	int numTriangles = getNumberOfTriangles(surfFile);
	int * Topology = getTopology(surfFile, numTriangles);
	double * colors = getColors(surfFile, numPoints);
	
	bool eliminateInteriorCavities = false;
//	bool eliminateInteriorCavities = true;
	SurfaceObject * result = new SurfaceObject(numPoints, VectorsAndNormals, numTriangles, Topology, eliminateInteriorCavities);
	
	result->addColors(colors);
	
	delete[](VectorsAndNormals);
	delete[](Topology);
	delete[](colors);

	fclose(surfFile);

	return result;
}


int getNumberOfPoints(FILE * surfaceFile)
{
	char * line = new char[1000];
	char * ptr;

	////get through the comments
	while( fgets(line, 1000, surfaceFile) != NULL){
		if(line[0] == '#'){ continue; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ continue; }
		if(strcmp(ptr, "GEOMETRY:") == 0){
			break;
		}
		else{
			printf("ERROR: Expected GEOMETRY: line that defines number of points in the surface.\n");
			printf("       Instead, found a line starting with this: [%s].  Exitting.\n", ptr);
			exit(1);
		}
	}

	////if the loop has ended, then we are at the GEOMETRY: line.
	////find out how many points there are.
	ptr = strtok(NULL, " \n\t\r");

	int result = atoi(ptr);

	delete[](line);

	return result;
}

double * getSurfaceGeometryAndNormals(FILE * surfaceFile, int numPts)
{
	int i = 0;
	char * line = new char[1000];
	char * ptr;
	double * result = new double[6*numPts];
	for(i = 0; i<numPts; i++){
		fgets(line, 1000, surfaceFile);

/*
		///if the line is a comment, warn and continue;
		bool iscomment = false;
		if(line[0] == '#'){ iscomment = true; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ iscomment = true; }
		if(iscomment){
			printf("WARNING: COMMENTS DETECTED INSIDE GEOMETRY: DEFINITION.  IGNORING\n");
			i = i - 1;
			continue;
		}
*/

		///otherwise its not a comment, so parse it.  Use linecopy.
		ptr = strtok(line, " \n\t\r");
		result[6*i+0] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+1] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+2] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+3] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+4] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[6*i+5] = atof(ptr);
	}

	delete[](line);

	return result;
}

int getNumberOfTriangles(FILE * surfaceFile)
{
	char * line = new char[1000];
	char * ptr;

	////get through the comments
	while( fgets(line, 1000, surfaceFile) != NULL){
		if(line[0] == '#'){ continue; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ continue; }
		if(strcmp(ptr, "TOPOLOGY:") == 0){
			break;
		}
		else{
			printf("ERROR: Expected GEOMETRY: line that defines number of points in the surface.\n");
			printf("       Instead, found a line starting with this: [%s].  Exitting.\n", ptr);
			exit(1);
		}
	}

	////if the loop has ended, then we are at the GEOMETRY: line.
	////find out how many points there are.
	ptr = strtok(NULL, " \n\t\r");

	int result = atoi(ptr);

	delete[](line);

	return result;
}

int * getTopology(FILE * surfaceFile, int numTris)
{
	int i = 0;
	char * line = new char[1000];
	char * ptr;
	int * result = new int[3*numTris];
	for(i = 0; i<numTris; i++){
		fgets(line, 1000, surfaceFile);

/*
		///if the line is a comment, warn and continue;
		bool iscomment = false;
		if(line[0] == '#'){ iscomment = true; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ iscomment = true; }
		if(iscomment){
			printf("WARNING: COMMENTS DETECTED INSIDE GEOMETRY: DEFINITION.  IGNORING\n");
			i = i - 1;
			continue;
		}
*/

		///otherwise its not a comment, so parse it.  Use linecopy.
		ptr = strtok(line, " \n\t\r");
		result[3*i+0] = atoi(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+1] = atoi(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+2] = atoi(ptr);
	}

	delete[](line);

	return result;	


}


///////////////////////////////////////////
///Gets the colors out of the SURF file.
double * getColors(FILE * surfaceFile, int numColors)
{
	double * result = new double[3*numColors];
	int i = 0;

	char * line = new char[1000];
	char * ptr;

	bool foundColors = false;

	////get through the comments
	while( fgets(line, 1000, surfaceFile) != NULL){
		if(line[0] == '#'){ continue; }
		ptr = strtok(line, " \n\t\r");
		if(ptr == NULL){ continue; }
		if(strcmp(ptr, "COLORS:") == 0){
			foundColors = true;
			break;
		}
		else{
			printf("ERROR: Expected COLORS: line that defines number of points in the surface.\n");
			printf("       Instead, found a line starting with this: [%s].  Exitting.\n", ptr);
			exit(1);
		}
	}
	
	if(!foundColors){
		delete[](result);
		delete[](line);
		return NULL;
	}
	

	///parse the colors
	for(i = 0; i<numColors; i++){
		fgets(line, 1000, surfaceFile);
		ptr = strtok(line, " \n\t\r");
		result[3*i+0] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+1] = atof(ptr);
		ptr = strtok(NULL, " \n\t\r");
		result[3*i+2] = atof(ptr);
	}
	
	delete[](line);
	return result;
}











