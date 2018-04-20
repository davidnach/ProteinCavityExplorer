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
 * File: StdAfx.cpp
 *       universal header file
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

#if !defined(_STDAFX_H_)
#define _STDAFX_H_

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////        INCLUDES        //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

/*--------------------- Libraries ---------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <limits.h>
#include <fstream>

/*------------------- Set Library ---------------------*/
#include "set.h"
#include "prime.h"
#include "defs.h"
#include "mathlib.h"
#include "funclib.h"

/*--------------- Atom Radius Lookup -------------------*/
#include "phimapGrid.h"

/*--------------- Martching Cubes headers------------------*/
#include "CubeTable.h"
#include "vertTriSet.h"

/*--------------- Surface File Input ------------------*/
#include "SurfaceObject.h"
#include "GeometryParser.h"

/*--------------- Lattice Hashing Headers ------------------*/
#include "LatticeHash.h"	//not related to targetGrid

/*------------------- Triangle Library ---------------------*/
#include "triIntersect.h"
#include "pointTriDist.h"

/*--------------- Surface Data Analysis-----------------*/
#include "LatticeObj.h"
#include "TriangleLatticeHash.h"

/*--------------- File gen and output -----------------*/
#include "SurfaceReconstruction.h"
#include "primativesProcessing.h"

/*--------------- Multiple Surface Processing-----------------*/
#include "SurfaceOutput.h"




////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////        DEFINITIONS        /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

/*------------Text Parsing #defs ----------*/

#define VERS_NUM 1.0

#define NUMCOLUMNS 10000					////parsing param for max # of columns
#define SMALL_NUM  0.000001 
#define NUM_INITIAL_RAY_TESTS 15
#define NUM_THOROUGH_RAY_TESTS 100

///Forces the cache to be large enough that we always only check once, so its never inconsistant.
///Forces the edge testing to generate a point even if we dont find one, so it snever inconsistant.
#define FORCE_LOGICAL_CONSISTANCY


#define HIGHLIGHT_INTERIOR 100
#define HIGHLIGHT_EXTERIOR 200
#define NOT_HIGHLIGHTED    000

#define SURFACE_BOOLEAN_UNION 100
#define SURFACE_BOOLEAN_INTERSECT 200
#define SURFACE_BOOLEAN_DIFFERENCE 300

#endif



