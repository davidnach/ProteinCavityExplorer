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
 * File: LatticeHash.cpp
 *       Spatial Hashing Structure interface
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

#ifndef LATTICE_HASH_H
#define LATTICE_HASH_H

#include "StdAfx.h"


// #####################################################################################################################################
// #####################################################################################################################################
// ####################    ###########    ######            ##            ##        #####       #####         ##########################
// ####################    ##########      #####            ##            ##        ###           ###         ##########################
// ####################    #########        ########    ##########    ########    ####    #####    ##    ###############################
// ####################    ########    ##    #######    ##########    ########    ####    ###########    ###############################
// ####################    #######    ####    ######    ##########    ########    ####    ###########       ############################
// ####################    #######    ####    ######    ##########    ########    ####    ###########       ############################
// ####################    #######            ######    ##########    ########    ####    ###########    ###############################
// ####################    #######            ######    ##########    ########    ####    ###########    ###############################
// ####################    #######    ####    ######    ##########    ########    ####    #####    ##    ###############################
// ####################         ##    ####    ######    ##########    ######        ###           ###         ##########################
// ####################         ##    ####    ######    ##########    ######        #####       #####         ##########################
// #####################################################################################################################################
// #####################################################################################################################################

// this code uses the following cube layout, with 3D axes to the right.
//       3---7            ^
//       |\  |\       \   |
//       | 1---5       +y |
//       | | | |        \ z+
//       2-|-6 |         \|  
//        \|  \|          0--x+-->    
//         0---4
// For sets specifying cube neighbors, they are specified [high-x][low-x][high-y][low-y][high-z][low-z]
// In the case where the cube has no neighbor in a certain position, NULL is provided in that space.
// Thus, sets specifying cube neighbors always have exactly six elements, but may have as few as three
// non-NULL members.

// Design Note:
// This class takes a 3D bounding box and a resolution, and creates a lattice of boxes that can be
// rapidly searched.
//

class primativesSet;

class LatticeHash
{
	public:
	
	// constructor and destructor
	LatticeHash( );
	virtual ~LatticeHash( );

	// setup for primativesSets;
	void setup( primativesSet * primSet, int resolution );

	// get objects in the same cube as a point, appends them to list
	set_t getNearbyObjs( double * pt, set_t list );

	// get objects in the same cubes as a segment, appends them to list
	set_t getNearbyObjs( double * v1, double * v2, set_t list );



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

	double res;							//resolution of the lattice
	double xneg, xpos, yneg, ypos, zneg, zpos;	//bounds of the lattice.
	int xdim, ydim, zdim;					//number of cubes on each side.
	
	set_t grid;							//This holds a list of single-dimensionally indexed pointers that represent
										//the boxes in the list.  They each point to a set object that consists of 
										//a list of indexes, which are the integer indexes of objects.  


};

#endif




