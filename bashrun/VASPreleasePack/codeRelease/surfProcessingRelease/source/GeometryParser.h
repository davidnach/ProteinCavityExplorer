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
 *       SURF file parser interface
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

#ifndef GEOMETRYPARSER_H
#define GEOMETRYPARSER_H

#include "StdAfx.h"

class SurfaceObject;

typedef struct Vertex {
  float x,y,z;             /* the usual 3-space position of a vertex */
} Vertex;

typedef struct Face {
  unsigned char intensity; /* this user attaches intensity to faces */
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
} Face;

///Parse the geometry file
SurfaceObject * parseGeometryFile(char * fileName); 

///find out how many points on the surface
int getNumberOfPoints(FILE * surfaceFile);

///get the vectors and surface normals on the surface
double * getSurfaceGeometryAndNormals(FILE * surfaceFile, int numPts);

///get the number of triangles
int getNumberOfTriangles(FILE * surfaceFile);

///get the triangles
int * getTopology(FILE * surfaceFile, int numTris);

///gets the colors
double * getColors(FILE * surfaceFile, int numColors);




#endif



