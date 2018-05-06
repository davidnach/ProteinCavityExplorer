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
 * Written by 
 *       Brian Y. Chen <chen@lehigh.edu>
 *
 * WWW URL: http://cse.lehigh.edu/~chen/
 * Email: chen@lehigh.edu
 * Documentation can be found here: 
 * http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
 *
 ***************************************************************************/

// SurfaceOutput.h : Header file for SurfaceOutput.cpp
//  which generates pdb surface files using the modified 
//  Troll library.

#ifndef _SURF_PROCESSING_H_
#define _SURF_PROCESSING_H_

#include "StdAfx.h"

class SurfaceObject;


// #################################################################################################################
// ###    ##  ##  #     ##     ###     #    #  ###     ###     ####   ###     ###    ##    #  ##  ##    ############
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  ##  ##     ##  ##  #  ##  ##  ##   #  #  ##  ###########
// ##  #####  ##  #  ##  #  ######  #####  ##  ###  ######  ##  ##  #  ##  ##  #  ######  ##      #  ###############
// ###  ####  ##  #    ###    ####    ###  ##  ###    ####  ##  #  ###  #    ####  #####  ##  #   #  ###############
// #####  ##  ##  #  #  ##  ######  #####  ##  ###  ######     ##       #  #  #####  ###  ##  ##  #  #   ###########
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  #####  ###  #  ##  #  ##  ##  ##  ##  #  ##  ###########
// ##  ##  #      #  ##  #  ######  #####  ##  ###  ######  #####  ###  #  ##  #  ##  ##  ##  ##  #  ##  ###########
// ###    ###    ##  ##  #  ######  ####    #    #     ###  #####  ###  #  ##  ##    ##    #  ##  ##    ############
// #################################################################################################################

// Top level parsing control.  Run this function to parse a given 
// file name into a SurfaceObject.  It calls the functions below.
SurfaceObject * parseGeometryFile(char * fileName); 

// find out how many points on the surface
int getNumberOfPoints(FILE * surfaceFile);

// get the vectors and surface normals on the surface
double * getSurfaceGeometryAndNormals(FILE * surfaceFile, int numPts);

// get the number of triangles
int getNumberOfTriangles(FILE * surfaceFile);

// get the triangles
int * getTopology(FILE * surfaceFile, int numTris);

// gets the colors
double * getColors(FILE * surfaceFile, int numColors);

// #################################################################################################################
// ###    ##  ##  #     ##     ###     #    #  ###     ####    ##  ##  #      #     ##  ##  #      #################
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  ##  #  ##  ###  ###################
// ##  #####  ##  #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  ##  #  ##  ###  ###################
// ###  ####  ##  #    ###    ####    ###  ##  ###    ####  ##  #  ##  ###  ###     ##  ##  ###  ###################
// #####  ##  ##  #  #  ##  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  #####  ##  ###  ###################
// ##  ##  #  ##  #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  #####  ##  ###  ###################
// ##  ##  #      #  ##  #  ######  #####  ##  ###  ######  ##  #  ##  ###  ###  #####  ##  ###  ###################
// ###    ###    ##  ##  #  ######  ####    #    #     ####    ###    ####  ###  ######    ####  ###################
// #################################################################################################################

// generateSURF:  This function generates a SURF file based on a surfaceObject and a set of
//                point indices that reference a patch on the surface.  The output is a
//                SURF file that contains only the points relevant to render the patch.
//                If the set of points is NULL, the entire surface is returned.
//                
//                =========================================================
//                  OUTPUT FORMAT:  A SURF file 
//                     -Anywhere in the file, empty lines and line starting with '#' are
//                     considered comments, and are ignored.
//                     -File has 2 sections: geometry and topology
//                
//                     The geometry section is always first, and starts with:
//                     GEOMETRY: <int>
//                     where <int> is the number of points to be specified.
//                     then there are <int> lines as follows:
//                     <float> <float> <float> <float> <float> <float>
//                     which stand for x,y,z, and xnormal, ynormal, znormal.
//                
//                     The topology section is always second, and starts with:
//                     TOPOLOGY: <int>
//                     where <int> specifies the number of triangles on the surface
//                     Then there are <int> lines as follows:
//                     <int> <int> <int>
//                     Each int stands for the 3 corners of the triangle
//                     and is an index into the array of points provided in the
//                     Geometry section.  I.e. the geometry section is indexed
//                     starting at zero, and ending at (size-1).
//                =========================================================
void generateSURF(SurfaceObject * surf, set_t points, char * outputFileName);






#endif



/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/




