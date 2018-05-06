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
 * File: phimapGrid.h
 *       parser for INSIGHT/DELPHI files
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

#ifndef PHIMAP_GRID_H
#define PHIMAP_GRID_H

#include "StdAfx.h"


class phimapGrid
{

public:
	phimapGrid( char * gridFileName );
	virtual ~phimapGrid();

	///gets the potential at a point
	double getPotential(double x, double y, double z);

	///variables
	float * midpoint;
	float * origin;
	float * extent;
	float id,ed,salt,perfil,scale,pr;
	float * pot;
	int bc;
	int npoints;	//total number of grid points (3 dims)
	int igrid;	//the number of grid points in one dimension
	float reaction_field_energy;
	int slicesize;
	float gridsize;	//physical size of one cube
	unsigned char state;


private:
	///reads the grid file
	void LoadInsight(char * fn);




};


#endif





