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
 * File: SurfaceOutput.h
 *       interface for outputting triangle meshes
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

#ifndef _SURFACE_OUTPUT_H_
#define _SURFACE_OUTPUT_H_

#include "StdAfx.h"


//////////////////////////////////////////////////////////
////Surface Generation based on values in INSIGHT grid
////
////This is located here because we need the trollbase stuff for it.  Otherwise it would be 
////located in simplifySurf
////
////Helper function Interpolates between two points with two charges, tells where the threshold point is.
double * pInterpolate(double c1, double * p1, double c2, double * p2, double thresh);
////Helper function returns true if inside, false if outside the desired region.
bool pCheck(double eVal, double threshold, int rangeMode);
////
////Surface Generation based on values in INSIGHT grid
////pdbfile, insight gridfile, threshold charge for surfface generation, resolution.
SurfaceObject * generatePhimapSurface(char * pdbFile, char * insightGrid, double thresh, int rangeMode, double res);
////
//////////////////////////////////////////////////////////


////generateSURF:  This function generates a SURF file based on a surfaceObject and a set of
////                  point indices that reference a patch on the surface.  The output is a
////                  SURF file that contains only the points relevant to render the patch.
////                  If the set of points is NULL, the entire surface is returned.
////                  
////                  =========================================================
////                    OUTPUT FORMAT:  A SURF file 
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
void generateSURF(SurfaceObject * surf, set_t points, char * outputFileName);









#endif




