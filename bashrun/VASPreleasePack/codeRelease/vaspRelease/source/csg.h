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
 * File: csg.h
 *       interface for CSG operations
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

#ifndef UNIFIED_MARCHING_CUBES_CSG
#define UNIFIED_MARCHING_CUBES_CSG

#define CSG_INTERSECTION 1000
#define CSG_DIFFERENCE   2000
#define CSG_UNION        3000

#include "StdAfx.h"

class Lattice;


// #####################################################################################################################################
// #####################################################################################################################################
// ######################       ########      ########       ###########################################################################
// ####################           ####          ####           #########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ###################    ###########    ##########    #################################################################################
// ###################    ############     ########    #################################################################################
// ###################    ##############      #####    ###      ########################################################################
// ###################    ##################    ###    ###      ########################################################################
// ###################    ###################    ##    #####    ########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ####################           ####          ####           #########################################################################
// ######################       ########      ########       ###########################################################################
// #####################################################################################################################################
// #####################################################################################################################################

class CSG
{
	public:
	CSG( SurfaceObject * surf1, SurfaceObject * surf2 );
	virtual ~CSG();
	
	// executive control function
	SurfaceObject * runCSG( int CSG_op, double resolution );


	// ##########################################################################################################################
	// ###     ###    ###     #########     ###     ##  ###  #      #     ###     ##  ###########################################
	// ##  ###  #  ##  #  ###  #######  ###  #  ###  #   ##  ###  ###  ##  #  ###  #  ###########################################
	// ##  ######  #####  ###  #######  ###  #  ###  #    #  ###  ###  ##  #  ###  #  ###########################################
	// ##  #######  ####  ############  ######  ###  #  #    ###  ###  ##  #  ###  #  ###########################################
	// ##  #########  ##  ############  ######  ###  #  ##   ###  ###    ###  ###  #  ###########################################
	// ##  ##########  #  ##   #######  ######  ###  #  ###  ###  ###  #  ##  ###  #  ###########################################
	// ##  ##########  #  ###  #######  ###  #  ###  #  ###  ###  ###  ##  #  ###  #  ###########################################
	// ##  ###  #  ##  #  ###  #######  ###  #  ###  #  ###  ###  ###  ##  #  ###  #  ###########################################
	// ###     ###    ###     #########     ###     ##  ###  ###  ###  ##  ##     ##     ########################################
	// ##########################################################################################################################
	
	// Functions for CSG capabilities.
	
	// Synchronize is a function that ensures that the lattices in both this lattice and
	// "lat", supplied as input, have the same dimensions.  This ensures that point and
	// cube indices one lattice map to the identical thing in the other lattice, so that
	// we can make the interior/exterior assessments we need to make in both without 
	// losing track of identity.  Call this function after setBounds() is called in both
	// lattices, NEVER before.  op is the CSG operation to be used.
	// This function also sets local variables (dimensions, resolution) for later use.
	void synchronize( Lattice * lat1, Lattice * lat2 );
	
	// This function generates the actual CSG surface.  It can only be called once sign
	// generation on two synchronized lattices has been computed.  After that is done
	// finalizeCSGsurface() finds all the nonempty cubes in the output surface (using
	// helper function csgNonemptyCubeTest()) and then generates surface on them in the 
	// same style as Lattice::finalizeSurface().  This data is then fed into the 
	// SurfaceObject constructor, and the SurfaceObject is generated and returned.
	SurfaceObject * finalizeCSGsurface( );
	
	// This function processes the case where either surface actually contains no geometry.
	// In these cases, we may return empty, or we may return one of the surfaces, based on
	// the operation and on the inputs.  
	SurfaceObject * degenerateProcessing( );
	
	
	
	// ####################################################################################################################################################
	// ###     ###    ###     #####  ###  #      #  #####     ##      #     #####      #  ##  #  ###  ##    ##      #    ##    ##  ###  ##    #############
	// ##  ###  #  ##  #  ###  ####  ###  #  #####  #####  ##  #  #####  ##  ####  #####  ##  #   ##  #  ##  ###  ####  ##  ##  #   ##  #  ##  ############
	// ##  ######  #####  ###  ####  ###  #  #####  #####  ##  #  #####  ##  ####  #####  ##  #    #  #  #######  ####  ##  ##  #    #  #  ################
	// ##  #######  ####  #########       #    ###  #####  ##  #    ###  ##  ####     ##  ##  #  #    #  #######  ####  ##  ##  #  #    ##  ###############
	// ##  #########  ##  #########  ###  #  #####  #####     ##  #####    ######  #####  ##  #  ##   #  #######  ####  ##  ##  #  ##   ####  #############
	// ##  ##########  #  ##   ####  ###  #  #####  #####  #####  #####  #  #####  #####  ##  #  ###  #  #######  ####  ##  ##  #  ###  #####  ############
	// ##  ##########  #  ###  ####  ###  #  #####  #####  #####  #####  ##  ####  #####  ##  #  ###  #  #######  ####  ##  ##  #  ###  #####  ############
	// ##  ###  #  ##  #  ###  ####  ###  #  #####  #####  #####  #####  ##  ####  #####  ##  #  ###  #  ##  ###  ####  ##  ##  #  ###  #  ##  ############
	// ###     ###    ###     #####  ###  #      #      #  #####      #  ##  ####  ######    ##  ###  ##    ####  ###    ##    ##  ###  ##    #############
	// ####################################################################################################################################################
	
	// Functions that support CSG capabilities
	
	// This function determins if a point is inside or outside based on the
	// desired CSG operation, and the interior/exterior state of the cubes
	// in each Lattice object.
	// This returns only 1 or 0, inside or outside.
	int isInsideCSG( int pointIndex );

	// This function determines if a given cube is nonempty in the output
	// surface, given the two input surfaces and the desired CSG operation.
	bool csgNonemptyCubeTest( int cubeIndex );

	// This function searches the nested set structure "edges" to determine if it contains
	// an index corresponding to c0 and c1 for a geometric point in edgePoints.  
	// If if a point exists, the function returns the index, otherwise it returns -1.
	int getEdgePointCSG(int p0, int p1);
	
	// This function inserts an edge into "edges", for indices c0 and c1, that corresponds
	// to point "pt".  "pt" itself is stored in edgePoints, using the index stored in edges.
	// Returns the index into edgePoints for this edge POint.
	int storeEdgePointCSG(int c0, int c1, double * pt );
	
	// This function computes the point at which the CSG boundary intersects the edge of
	// the cube.  This is achieved by first getting the interior/exterior status of the
	// points in each Lattice, then making a decision based on the appropriate intersection
	// point.
	double * computeEdgeCrossingPointCSG( int p0, int p1 );



	// ##########################################################################################################################
	// ##  ###  #### ####     ##    #### ####     ##  #####      ##    ##########################################################
	// ##  ###  ###   ###  ##  ##  ####   ###  ##  #  #####  #####  ##  #########################################################
	// ##  ###  ##     ##  ##  ##  ###     ##  ##  #  #####  #####  #############################################################
	// ##  ###  #   #   #  ##  ##  ##   #   #     ##  #####    ####  ############################################################
	// ##  ###  #  ###  #    ####  ##  ###  #  ##  #  #####  ########  ##########################################################
	// ##  ###  #       #  #  ###  ##       #  ##  #  #####  #########  #########################################################
	// ###  #  ##  ###  #  ##  ##  ##  ###  #  ##  #  #####  #########  #########################################################
	// ####   ###  ###  #  ##  ##  ##  ###  #  ##  #  #####  #####  ##  #########################################################
	// ##### ####  ###  #  ##  #    #  ###  #     ##      #      ##    ##########################################################
	// ##########################################################################################################################

	//Variables/////////////////////////////////////////////////////////////////

	int op;								// The CSG operation that we will execute.

	bool runnable;							// constructor sets this to false if either input 
										// surface is empty.  True otherwise.

	Lattice * lat1;						// Lattices that we will use for the CSG operation.
	Lattice * lat2;

	double res;							// resolution of the lattice
	double xneg, xpos, yneg, ypos, zneg, zpos;	// bounds of the lattice.
	int xdim, ydim, zdim;					// number of cubes on each side.

	set_t edges;							// This set contains sets that map to a point in space.
	set_t edgePoints;						// The top level set is indexed by points, and the
										// lower level set is also indexed by points, so you get
										// a way to index point-to-point connections - i.e. edges.
										// The lower level set contains int* values (int[1]) that index
										// into edgePOints, where you can get the actual double vector


};

// #####################################################################################################################################
// #####################################################################################################################################
// ######################       ########      ########       ###########################################################################
// ####################           ####          ####           #########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ###################    ###########    ##########    #################################################################################
// ###################    ############     ########    #################################################################################
// ###################    ##############      #####    ###      ########################################################################
// ###################    ##################    ###    ###      ########################################################################
// ###################    ###################    ##    #####    ########################################################################
// ###################    #####    ##    ####    ##    #####    ########################################################################
// ####################           ####          ####           #########################################################################
// ######################       ########      ########       ###########################################################################
// #####################################################################################################################################
// #####################################################################################################################################


#endif




/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/


