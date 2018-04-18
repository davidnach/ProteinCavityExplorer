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
 * File: StdAfx.h
 *       include file for standard system include files,
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
#include <time.h>
#include <limits.h>
#include <math.h>


/*------------------- Set Library ---------------------*/
#include "set.h"
#include "prime.h"
#include "defs.h"

/*--------------- libraries by Brian Chen ------------------*/
#include "mathlib.h"
#include "triIntersect.h"
#include "progressBar.h"

/*--------------- Martching Cubes case-logic class ------------------*/
#include "CubeTable.h"

/*--------------- Data Representation and Maintenance Headers ------------------*/
#include "SurfaceObject.h"

/*--------------- File Handling Headers ------------------*/
#include "surfProcessing.h"


/*--------------- Marching Cubes Simplification and Data Structure ------------------*/
#include "Lattice.h"
#include "csg.h"


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////        DEFINITIONS        /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

#define SMALL_NUM  0.000001 				////an epsilon value for some math functions
#define VERS_NUM   1.0

#endif





/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/





