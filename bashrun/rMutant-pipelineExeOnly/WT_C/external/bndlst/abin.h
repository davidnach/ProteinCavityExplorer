/* name: abin.c                                                 */
/* author: J. Michael Word (port from dcr and mez fortran code) */
/* date written: 2/20/96                                        */
/* purpose: organize atoms into bins to aid finding neighbors   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef ABIN_H
#define ABIN_H 1

#include <math.h>
#include "geom3d.h"

typedef struct {  /* 3D bounding region */
	point3d min, max;
} region;

#define ATOMS_IN_THE_SAME_RES(aa, bb)         \
 ( ((aa)->rescnt == (bb)->rescnt)             \
   || (  (aa)->resid      == (bb)->resid      \
      && (aa)->resInsCode == (bb)->resInsCode \
      && (aa)->chain      == (bb)->chain      \
      && (aa)->model      == (bb)->model) )

#define FPRINT_ATOMNAME(outf, aptr) \
              fprintf((outf), "%4.4s%c%3.3s %c%4d%c",            \
                           (aptr)->atomname, (aptr)->altConf,    \
                           (aptr)->resname, (aptr)->chain,       \
                           (aptr)->resid, (aptr)->resInsCode)

typedef struct atom_t {
	struct atom_t *next;      /* primary list of atoms         */
	struct atom_t *nextInBin; /* list of atoms in same bin     */
	struct atom_t *scratch;   /* used for arbitrary atom lists */

	char* bondedto; /* for H atoms where the parent atom is specified  */
	                /* for heavy atoms, the list of bonded heavy atoms */
			/* are listed (for standard residues)              */

	int mark; /* general purpose (e.g. used for marking bonded atoms) */

	long atomCount;  /* id number */

        int  flags;      /* selection mechanism */
	int  props;      /* property flags      */

	int  file;       /* which file did atom come from? */
	int  model;      /* which model is the atom for? */
	int  resid;      /* residue number */
	int  rescnt;     /* residue order 1, 2, 3, ... */

	int  elem;       /* elememt number  */
	int  atomclass;  /* either elem or (for H) the parent elem or -1 (for unknown) */

	point3d loc;     /* xyz position    */
	int ix, iy, iz;  /* address in bins */

	float radius;    /* VDW radius      */
	float covRad;    /* covalent radius */

	float occ;       /* fractional occupancy */
	float bval;      /* temperature factor   */

	char  segid[5];  /* segment identifier   */

	char    resname[5];
	char   atomname[5];

	char chain;      /* peptide chain code */
	char resInsCode; /* insertion code */
	char altConf;    /* atom alternate conformation code */
} atom, *atomPtr;

typedef struct {
	point3d min, max, sz; /* bounding box and span */
	int nx, ny, nz;       /* bin size in each dimension */
	float  delta;         /* size bin edges in angstroms */
	atomPtr ***list;      /* 3d matrix of atom lists */
} atomBins;

#define OFFSET(val, min, delta) floor(1.0 + ((val)-(min))/(delta))

atomBins* initBins(region *boundingBox, float delta);
void disposeBins(atomBins* bins);
void updateBoundingBox(point3d *loc, region *boundingBox);
void addNeighbor(atom *a, atomBins *bins);
void binLoc(point3d *loc, int *ix, int *iy, int *iz, atomBins *bins);

#endif

