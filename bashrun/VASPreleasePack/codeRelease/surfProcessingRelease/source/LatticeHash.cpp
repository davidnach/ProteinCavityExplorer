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
 *       Spatial Hashing Structure implementation
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

#include "LatticeHash.h"


/////////////////////////////////////////////////////////////////////////////////////
// Nothing is computed in the constructor.  Just set up the data.
/////////////////////////////////////////////////////////////////////////////////////
LatticeHash::LatticeHash(  )
{
	res = 0;							// null resolution
	
	///set the initial boundary values
	xneg = HUGE_VAL;	xpos = -HUGE_VAL;
	yneg = HUGE_VAL;	ypos = -HUGE_VAL;
	zneg = HUGE_VAL;	zpos = -HUGE_VAL;

	///set the initial grid extents.
	xdim = 0;
	ydim = 0;
	zdim = 0;
	
	///setup the grid
	grid = alloc_set(SP_MAP);	

}


/////////////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////////////
LatticeHash::~LatticeHash()
{
	int i = 0;
	for(i = 0; i<size_set(grid); i++){
		set_t temp = (set_t) mapsto_set(grid, grid[i]); 
		free_set(temp);
	}
	
	free_set(grid);
	
}




/////////////////////////////////////////////////////////////////////////////////////
// Set up a lattice hash with a primSet.
/////////////////////////////////////////////////////////////////////////////////////
void LatticeHash::setup( primativesSet * primSet, int resolution ){
	
	//set the resolution
	res = resolution;

	//get the bounds	
	double * tempBounds = primSet->getPrimSetBounds();
	
	///read out the bounds from the primSet
	xneg = tempBounds[0];
	xpos = tempBounds[1];
	yneg = tempBounds[2];
	ypos = tempBounds[3];
	zneg = tempBounds[4];
	zpos = tempBounds[5];
	delete[](tempBounds);

	///random nudge to offset axis aligned geometry.
	double xnudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double ynudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	double znudge = ( ( (double) (random()%1000) ) / 1000.0 ) * res;
	xneg -= xnudge;
	yneg -= ynudge;
	zneg -= znudge;

	////////////////////////////////////////
	printf("LATTICEHASH DIMENSIONS:            xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	////////////////////////////////////////

	//padding on low side.
	xneg -= res;	yneg -= res;	zneg -= res;

	///+1 for fractional cube round up (int conversion truncates)
	///+1 for high side padding
	///= +2 total additions.
	xdim = ((int) ((xpos-xneg)/res)) + 2; 	///this is the total number of CUBES.  Not lines.
	ydim = ((int) ((ypos-yneg)/res)) + 2;
	zdim = ((int) ((zpos-zneg)/res)) + 2;
	///Reset the high dimensions based on the sizes of the cubes, since there is fractional roundup.
	xpos = xneg + xdim*res;
	ypos = yneg + ydim*res;
	zpos = zneg + zdim*res;
	
	////////////////////////////////////////
	printf("EXPANDED LATTICEHASH DIMENSIONS:   xneg: %f xpos: %f yneg: %f ypos: %f zneg: %f zpos: %f\n", xneg, xpos, yneg, ypos, zneg, zpos);
	printf("LATTICEHASH CUBE DIMENSIONS:       xdim: %i, ydim: %i, zdim: %i\n", xdim, ydim, zdim);
	printf("now adding primitives from primativesSet...\n");
	////////////////////////////////////////

	int l = 0;
	for(l = 0; l<primSet->size(); l++){
//		printf("inserting\n");
		Primative * p = primSet->getPrim(l);
		double * tempBounds = p->getBounds();

		//find out which cubes the bounding box touches
		double sxmin = tempBounds[0];
		double sxmax = tempBounds[1];
		double symin = tempBounds[2];
		double symax = tempBounds[3];
		double szmin = tempBounds[4];
		double szmax = tempBounds[5];
		delete[](tempBounds);

		//get the index boundaries of the bounding box, as if the grid was infinite.
		int ihi, ilo, jhi, jlo, khi, klo;
		if( fmod( (sxmin-xneg), res ) == 0 ){ ilo = ((int)((sxmin-xneg)/res)) -1; }
		else{ ilo = ((int)((sxmin-xneg)/res)); }
		if( fmod( (symin-yneg), res ) == 0 ){ jlo = ((int)((symin-yneg)/res)) -1; }
		else{ jlo = ((int)((symin-yneg)/res)); }
		if( fmod( (szmin-zneg), res ) == 0 ){ klo = ((int)((szmin-zneg)/res)) -1; }
		else{ klo = ((int)((szmin-zneg)/res)); }
		ihi = (int)( (sxmax-xneg)/res );	// if it is on the exact high end of a cube, it will 
		jhi = (int)( (symax-yneg)/res );	//  give the whole additional cube, so the checks above
		khi = (int)( (szmax-zneg)/res );	//  do not apply to ihi, jhi, khi.

		// Check index bounds.  Breaking index bounds should NOT be possible
		// Unless the lattice was not built for this surface.
		if( ihi>=xdim || ilo<0 || jhi>=ydim || jlo<0 || khi>=zdim || klo<0 ){
			printf("ERROR [LatticeHash Bounds]: Attempted triangle insertion violates LatticeHash bounds.\n");
			printf("ERROR [LatticeHash Bounds]: This LatticeHash was not built using your primitives.  Exitting.\n");
			exit(1);
		}

		// retrieve the cubes.  Note that the cubes are indexed in three dimensions, so we must
		// first span the indices, then transform them into a single dimension index (see myindex)
		int i,j,k;
		for(i = ilo; i<=ihi; i++){
			for(j = jlo; j<=jhi; j++){
				for(k = klo; k<=khi; k++){
					// get the cube index
					int myIndex = (i*ydim*zdim) + (j*zdim) + k;
					
					// insert the index of this primivitve into the grid.
					set_t myGridBox;
					if( contains_set(grid, myIndex) ){ myGridBox = (set_t) mapsto_set(grid, myIndex); }
					else{
						myGridBox = alloc_set(0);
					}
					
					// store the index of the primitive (l) in the cube.
					myGridBox = put_set( myGridBox, l );
					
					// now put the box back into the grid.
					grid = associate_set( grid, myIndex, myGridBox );
				}
			}
		}
	}
	
	
}






/////////////////////////////////////////////////////////////////////////////////////
// get objects in the same cube as a point
/////////////////////////////////////////////////////////////////////////////////////
set_t LatticeHash::getNearbyObjs( double * pt, set_t list )
{
	//declare variables
	int ibox, jbox, kbox;
	double x = pt[0];
	double y = pt[1];
	double z = pt[2];

	//find the box that contains the point	
	if( fmod( (x-xneg), res ) == 0 ){ ibox = ((int)((x-xneg)/res)) -1; }
	else{ ibox = ((int)((x-xneg)/res)); }
	if( fmod( (y-yneg), res ) == 0 ){ jbox = ((int)((y-yneg)/res)) -1; }
	else{ jbox = ((int)((y-yneg)/res)); }
	if( fmod( (z-zneg), res ) == 0 ){ kbox = ((int)((z-zneg)/res)) -1; }
	else{ kbox = ((int)((z-zneg)/res)); }

	// Check index bounds.  Breaking index bounds should NOT be possible
	// Unless the lattice was not built for this surface.
	if( ibox>=xdim || ibox<0 || jbox>=ydim || jbox<0 || kbox>=zdim || kbox<0 ){
		
//		printf("OMG: [%f %f %f], %f-%f %f-%f %f-%f\n", x, y, z, xneg, xpos, yneg, ypos, zneg, zpos);
//		printf("WTF: [%i %i %i], %i-%i %i-%i %i-%i\n", ibox, jbox, kbox, 0, xdim, 0, ydim, 0, zdim);
		return list;
	}

	//get the array index of the box
	int myIndex = (ibox*ydim*zdim) + (jbox*zdim) + kbox;

	// get the box from the grid.
	set_t myGridBox;
	
	// if the grid contains this box (i.e. something nonzero here, then
	// insert the insides.  otherwise do nothing (i.e. no else)
	if( contains_set(grid, myIndex) ){ 
		myGridBox = (set_t) mapsto_set(grid, myIndex); 

		int i = 0;
		//go through the set and insert the related objects into the results;
		for(i = 0; i<size_set(myGridBox); i++){
			list = put_set(list, myGridBox[i]);
		}
	}

	//return the results in the list;
	return list;	
}



/////////////////////////////////////////////////////////////////////////////////////
// get objects in the same cubes as a segment, appends them to list
// based on getSegmentCubes in Lattice.cpp in VASP.
/////////////////////////////////////////////////////////////////////////////////////
set_t LatticeHash::getNearbyObjs( double * v1, double * v2, set_t list )
{
//	int i = 0;
//	int j = 0;
	
	set_t listOfCubes = alloc_set(0);
	int xind, yind, zind;

//	printf("RUNNING GETSEGMENTCUBES--------------------------------------------------------\n");
//	printf("xneg: %f  yneg: %f  zneg: %f   res: %f\n", xneg, yneg, zneg, res);
	
	double fmodx = fabs(fmod( v1[0]-xneg, res ));
	double fmody = fabs(fmod( v1[1]-yneg, res ));
	double fmodz = fabs(fmod( v1[2]-zneg, res ));

	if( fmodx<.00001 || fmody<.00001 || fmodz<.00001 ){
		//the point is on a corner or face or edge, so nudge closer to v1 so that it's inside the cube that the
		//segment will pass through, and find out what segment it is.
		double * norm = normalizeVector( v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2] );
		double nx = v1[0] + ((res*.01)*norm[0]);
		double ny = v1[1] + ((res*.01)*norm[1]);
		double nz = v1[2] + ((res*.01)*norm[2]);

		xind = (int) ((nx-xneg)/res);	//if(xind == xdim){ xind = xdim-1; }//these not necessary
		yind = (int) ((ny-yneg)/res);	//if(yind == ydim){ yind = ydim-1; }//since we are inside
		zind = (int) ((nz-zneg)/res);	//if(zind == zdim){ zind = zdim-1; }//the cube via nudge (nx,ny,nz)
		delete[](norm);
	}
	else{
//		printf("I'm not on a corner or face or edge\n");
		xind = (int) ((v1[0]-xneg)/res);	//if(xind == xdim){ xind = xdim-1; }//these not necessary
		yind = (int) ((v1[1]-yneg)/res);	//if(yind == ydim){ yind = ydim-1; }//because we dont
		zind = (int) ((v1[2]-zneg)/res);	//if(zind == zdim){ zind = zdim-1; }//need a nudge anyway
	}

	///put the first point in.
	listOfCubes = put_set(listOfCubes, xind*(ydim*zdim) + yind*zdim + zind );

	while( 1 ){
		int xbot=0, xtop=0, ybot=0, ytop=0, zbot=0, ztop=0;

		// declare the cube constraints
		double xlo = xneg + ((xind+0)*res);	double xhi = xneg + ((xind+1)*res);
		double ylo = yneg + ((yind+0)*res);	double yhi = yneg + ((yind+1)*res);
		double zlo = zneg + ((zind+0)*res);	double zhi = zneg + ((zind+1)*res);

		// find out what walls it hits, unless we are at a border cube.
		if( xind!=0      ){ xbot = segAAquadCollision(v1, v2, xlo, xlo, ylo, yhi, zlo, zhi, NULL); }
		if( xind!=xdim-1 ){ xtop = segAAquadCollision(v1, v2, xhi, xhi, ylo, yhi, zlo, zhi, NULL); }
		if( yind!=0      ){ ybot = segAAquadCollision(v1, v2, xlo, xhi, ylo, ylo, zlo, zhi, NULL); }
		if( yind!=ydim-1 ){ ytop = segAAquadCollision(v1, v2, xlo, xhi, yhi, yhi, zlo, zhi, NULL); }
		if( zind!=0      ){ zbot = segAAquadCollision(v1, v2, xlo, xhi, ylo, yhi, zlo, zlo, NULL); }
		if( zind!=zdim-1 ){ ztop = segAAquadCollision(v1, v2, xlo, xhi, ylo, yhi, zhi, zhi, NULL); }

		// eliminate cases where v1 is on a corner.
		if( xbot==1 && fabs(v1[0]-xlo)<.00001 ){ xbot = 0; }
		if( xtop==1 && fabs(v1[0]-xhi)<.00001 ){ xtop = 0; }
		if( ybot==1 && fabs(v1[1]-ylo)<.00001 ){ ybot = 0; }
		if( ytop==1 && fabs(v1[1]-yhi)<.00001 ){ ytop = 0; }
		if( zbot==1 && fabs(v1[2]-zlo)<.00001 ){ zbot = 0; }
		if( ztop==1 && fabs(v1[2]-zhi)<.00001 ){ ztop = 0; }

		//the segment hits the cube in two places, one of which comes from the cube we were just in.
		//that cube is already in the listOfCubes list - so check for it, and then remove it.
		if( xbot==1 && contains_set(listOfCubes, (xind-1)*(ydim*zdim) + (yind+0)*zdim + (zind+0) ) ){ xbot = 0; }
		if( xtop==1 && contains_set(listOfCubes, (xind+1)*(ydim*zdim) + (yind+0)*zdim + (zind+0) ) ){ xtop = 0; }
		if( ybot==1 && contains_set(listOfCubes, (xind+0)*(ydim*zdim) + (yind-1)*zdim + (zind+0) ) ){ ybot = 0; }
		if( ytop==1 && contains_set(listOfCubes, (xind+0)*(ydim*zdim) + (yind+1)*zdim + (zind+0) ) ){ ytop = 0; }
		if( zbot==1 && contains_set(listOfCubes, (xind+0)*(ydim*zdim) + (yind+0)*zdim + (zind-1) ) ){ zbot = 0; }
		if( ztop==1 && contains_set(listOfCubes, (xind+0)*(ydim*zdim) + (yind+0)*zdim + (zind+1) ) ){ ztop = 0; }

		//error check
		if( xbot+xtop+ybot+ytop+zbot+ztop>2 ){ 
			printf("ERROR [getSegmentCubes]: segment intersected %i cubes! WTF\n", xbot+xtop+ybot+ytop+zbot+ztop);
			if(xbot){ printf("xbot\n"); }
			if(xtop){ printf("xtop\n"); }
			if(ybot){ printf("ybot\n"); }
			if(ytop){ printf("ytop\n"); }
			if(zbot){ printf("zbot\n"); }
			if(ztop){ printf("ztop\n"); }
			
			printf("v1: %f %f %f\n", v1[0], v1[1], v1[2] );
			printf("v2: %f %f %f\n", v2[0], v2[1], v2[2] );
			printf("xlo %f xhi %f ylo %f yhi %f zlo %f zhi %f \n", xlo, xhi, ylo, yhi, zlo, zhi );
			//printf("counter = %i\n");
			
			//if(v1Corner){ printf("V1 CORNER true\n"); }
		}

		//if one of these is 1, then we can go in that dir bc it's not on the wall already.
		//becayse we checked that above)
		int listOfCubesSize = size_set(listOfCubes);
		if( xbot==1 ){ xind = xind-1; listOfCubes = put_set( listOfCubes, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( xtop==1 ){ xind = xind+1; listOfCubes = put_set( listOfCubes, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( ybot==1 ){ yind = yind-1; listOfCubes = put_set( listOfCubes, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( ytop==1 ){ yind = yind+1; listOfCubes = put_set( listOfCubes, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( zbot==1 ){ zind = zind-1; listOfCubes = put_set( listOfCubes, xind*(ydim*zdim) + yind*zdim + zind ); }
		if( ztop==1 ){ zind = zind+1; listOfCubes = put_set( listOfCubes, xind*(ydim*zdim) + yind*zdim + zind ); }

		if( size_set(listOfCubes) == listOfCubesSize ){ 
			break;
		}
		
	}
	
	///now we have the list of cubes.  Now get the objects out of the cubes.
	int i = 0;
	for(i = 0; i< size_set(listOfCubes); i++){

		//get next cube index from listOfCubes, and use that index 
		//to get the resulting list from grid.
		int myIndex = listOfCubes[i];
		if( contains_set(grid, myIndex) ){
			set_t myGridBox = (set_t) mapsto_set(grid, myIndex); 
	
			int j = 0;
			//go through the set and insert the related objects into the results;
			for(j = 0; j<size_set(myGridBox); j++){
				list = put_set(list, myGridBox[j]);
			}
		}


	
	}
	
	free_set(listOfCubes);
	
	//now return the final list.
	return list;

}














