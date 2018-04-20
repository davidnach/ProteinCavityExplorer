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
 * File: phimapGrid.cpp
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



#include "phimapGrid.h"





//phimapGrid::phimapGrid( char * pdbFileName, char * gridFileName )
phimapGrid::phimapGrid( char * gridFileName )
{

	//set the grid initially zero.
	pot=NULL;
	
	//allocate spatial bounds.
	origin = new float[3];
	origin[0] = 0;	origin[1] = 0;	origin[2] = 0; 
	extent = new float[3];
	extent[0] = 0;	extent[1] = 0;	extent[2] = 0; 
	midpoint = new float[3];
	midpoint[0] = 0;	midpoint[1] = 0;	midpoint[2] = 0; 

	///parse the centroid	
//	parsePdbForCentroid( pdbFileName );
	
	LoadInsight(gridFileName);
}

phimapGrid::~phimapGrid()
{
	delete[](origin);
	delete[](extent);
	delete[](midpoint);
	delete[](pot);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void phimapGrid::LoadInsight(char *fn)

{
	bc=0;
	salt=0;
	id=0;
	ed=0;
	perfil=0;
	pr=0;
//	atom=structure->GetSelection(allatoms);	?????

	///open the binary file handle, read only.
	int fh;
	fh=open(fn,O_RDONLY);

	//temp vars
	char * temp = new char[132];
	int i;
	float ext;
	float t,x0,x1,y0,y1,z0,z1;

//*************************************************************/
//*
//*  A note on Fortran 'write' records.
//*  Four bytes are always reserved at either end of a 
//*  record every time Fortran calls the function to 'write' 
//*  records to a file, as follows:
//*
//*    XXXX data1 data2 ... dataN XXXX 
//* 
//*  where each X represents a byte.  See 
//*  http://astronomy.swin.edu.au/~pbourke/dataformats/fortran/
//*  for more information.
//*
//*************************************************************/

	// qnifft write: toplbl
	read(fh,temp,60+8);        // toplbl + rec header/trailer

	// qnifft write: ivary,nbyte,intdat,extent*3,xang,yang,zang,
	//               xstart,xend,ystart,yend,zstart,zend,intx,inty,intz
	read(fh, temp, 4 );         // rec header
	read(fh,(char *)&i,4);     // ivary = 0
	read(fh,(char *)&i,4);     // nbyte = 4

	///error check.
	if (i!=4){
		printf("ERROR: cannot read phimap binary file\n");
		exit(1);
	} 

	read(fh,temp,4);           // intdat = 0
	read(fh,(char *)&ext,4);   // ext * 3
	read(fh,(char *)&ext,4);
	read(fh,(char *)&ext,4);
	read(fh,(char *)&t,4);     // angle (x,y,z) = 90,90,90
	read(fh,(char *)&t,4);
	read(fh,(char *)&t,4);

	read(fh,(char *)&x0,sizeof(float));  // x/y/z start/end
	read(fh,(char *)&x1,sizeof(float));
	read(fh,(char *)&y0,sizeof(float));
	read(fh,(char *)&y1,sizeof(float));
	read(fh,(char *)&z0,sizeof(float));
	read(fh,(char *)&z1,sizeof(float));

	///we are assigning igrid for the first time in this object
	read(fh,(char *)&igrid,sizeof(int)); // igrid * 3
	read(fh,(char *)&igrid,sizeof(int));
	read(fh,(char *)&igrid,sizeof(int)); 
	read(fh,temp, 4 );         // rec trailer 

	++igrid;

	// derivation of equations below 
	// follow from qnifft22 source code
	// (see also delphi documentation)

	float range = (x1-x0) * ext / 2.f;   

	//midpoint.x = x0 * ext + range; 
	midpoint[0] = x0 * ext + range; 
	//midpoint.y = y0 * ext + range; 
	midpoint[1] = y0 * ext + range; 
	//midpoint.z = z0 * ext + range; 
	midpoint[2] = z0 * ext + range; 

	//set scale
	scale = float (igrid-1) / (2.f*range); 

	//set reaction field energy to zero (set for the first time.
	reaction_field_energy = 0.f;

	//clear the memory in case. //this never gets called.  original code clears it earlier.
	if(pot!=NULL){
		delete pot;
	}

	// read fileblock into memory then trim
	// records using memmove; realloc to free
	// unused memory after

	//set npoints.
	npoints=igrid*igrid*igrid;
	int fpoints=igrid*igrid*(igrid+2);
	
	//allocate the grid
	pot=new float[fpoints];
	//read the grid
	read(fh,(char *)pot,fpoints*sizeof(float)); 

	///close the file.
	close(fh);

	// each igrid block is a fortran write record;
	// move the bytes to get rid of the headers and 
	// trailers of each record

	//now we allocate sliceSize
	slicesize=igrid*igrid;
	for (i=0; i<slicesize; ++i) {
		int targ = sizeof(float)*((igrid+2)*i+1);
		int dest = sizeof(float)*(igrid*i);
		memmove((char *)pot+dest, (char *)pot+targ, igrid*sizeof(float));
	}
	pot=(float *)realloc(pot,npoints*sizeof(float));

	// calc origin & extent

	///now we are allocating gridSize (physical size of one cube)
	gridsize=1.f/scale;

	///this is a vector	
//	Vector unit(gridsize,gridsize,gridsize);

	//this is a scalar
	//igrid is the number of grid points in one dimneiosn
	int midgrid=(igrid-1)/2;

	//origin=midpoint-unit*midgrid;
	
	origin[0] = midpoint[0] - (gridsize*midgrid);
	origin[1] = midpoint[1] - (gridsize*midgrid);
	origin[2] = midpoint[2] - (gridsize*midgrid);
	
	//extent=origin+unit*igrid;
	
	extent[0] = origin[0] + (gridsize*igrid);
	extent[1] = origin[1] + (gridsize*igrid);
	extent[2] = origin[2] + (gridsize*igrid);
	
	//now we are allocating state.
	state=1;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// Added by brian chen 04/12/2007
///// COPIES EXPLICITLY FROM operator[] - could not get that to work

double phimapGrid::getPotential(double x, double y, double z)
{
	int nx,ny,nz,nx1,ny1,nz1;
	float xgr,ygr,zgr,dx,dy,dz;
	float a1,a2,a3,a4,a5,a6,a7,a8;
	float phi;
	
	// Error if outside grid.
	if(x>extent[0]) {
		printf("TROLLBASE ERROR - OUT OF BOUNDS: Point (%f, %f, %f) is outside +x extent %f\n", x,y,z, extent[0]);
		exit(1);	//	throw(Error(TE_PhimapError,"point outside grid, potential undefined"));
	}
	if(y>extent[1]) {
		printf("TROLLBASE ERROR - OUT OF BOUNDS: Point (%f, %f, %f) is outside +y extent %f\n", x,y,z, extent[1]);
		exit(1);	//	throw(Error(TE_PhimapError,"point outside grid, potential undefined"));
	}
	if(z>extent[2]) {
		printf("TROLLBASE ERROR - OUT OF BOUNDS: Point (%f, %f, %f) is outside +z extent %f\n", x,y,z, extent[2]);
		exit(1);	//	throw(Error(TE_PhimapError,"point outside grid, potential undefined"));
	}
	if(x<origin[0]) {
		printf("TROLLBASE ERROR - OUT OF BOUNDS: Point (%f, %f, %f) is outside -x origin %f\n", x,y,z, origin[0]);
		exit(1);	//	throw(Error(TE_PhimapError,"point outside grid, potential undefined"));
	}
	if(y<origin[1]) {
		printf("TROLLBASE ERROR - OUT OF BOUNDS: Point (%f, %f, %f) is outside -y origin %f\n", x,y,z, origin[1]);
		exit(1);	//	throw(Error(TE_PhimapError,"point outside grid, potential undefined"));
	}
	if(z<origin[2]) {
		printf("TROLLBASE ERROR - OUT OF BOUNDS: Point (%f, %f, %f) is outside -z origin %f\n", x,y,z, origin[2]);
		exit(1);	//	throw(Error(TE_PhimapError,"point outside grid, potential undefined"));
	}
	
	// Calculate lower left bottom grid point
	dx=x-origin[0];
	dy=y-origin[1];
	dz=z-origin[2];
	
	// !!! casting next three floats to int per warning
	nx=(int)(dx*scale);
	nx1=nx+1;
	if(nx1>igrid) nx1=nx;
	
	ny=(int)(dy*scale);
	ny1=ny+1;
	if(ny1>igrid) ny1=ny;
	
	nz=(int)(dz*scale);
	nz1=nz+1;
	if(nz1>igrid) nz1=nz;
	
	xgr = dx-nx*gridsize;
	ygr = dy-ny*gridsize;
	zgr = dz-nz*gridsize;
	
	// Calculate coefficients of trilinear function.
	a8=pot[nz*slicesize+ny*igrid+nx];
	a7=pot[nz1*slicesize+ny*igrid+nx]-a8;
	a6=pot[nz*slicesize+ny1*igrid+nx]-a8;
	a5=pot[nz*slicesize+ny*igrid+nx1]-a8;
	a4=pot[nz1*slicesize+ny1*igrid+nx]-a8-a7-a6;
	a3=pot[nz1*slicesize+ny*igrid+nx1]-a8-a7-a5;
	a2=pot[nz*slicesize+ny1*igrid+nx1]-a8-a6-a5;
	a1=pot[nz1*slicesize+ny1*igrid+nx1]-a8-a7-a6-a5-a4-a3-a2;
	
	// Determine value of phi
	phi=a1*xgr*ygr*zgr+a2*xgr*ygr+a3*xgr*zgr+a4*ygr*zgr+a5*xgr+a6*ygr+a7*zgr+a8;
	
	return (double)phi;

}

// Added by brian chen 04/12/2007
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
  
  
  
  