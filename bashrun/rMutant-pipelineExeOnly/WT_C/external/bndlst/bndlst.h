/* name: bndlst.h                        */
/* author: J. Michael Word               */
/* date written: 4/9/99                  */
/* purpose: definitions for connection sw*/

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef BNDLST_H
#define BNDLST_H 1

#include <stdio.h>
#include "utility.h"
#include "abin.h"
#include "readPDBrecs.h"
#include "select.h"
#include "atomprops.h"
#include "stdconntable.h"

void bndlsthelp();
atom* processCommandline(int argc, char **argv, region *bbox,
   char **idString, int *maxBonds, int *noBonds,
   int *modelNum, float *probeRad);
void reverseAtomList(atom **allAtoms);
atom* loadAtoms(FILE *fp, atom *atomlist, region *bbox, int modelNum);
atomBins* binAtoms(atom *allAtoms, region *boundingBox, float probeRad);
float getRadius(int at);
atom * newAtom(char *rec, int file, int model);
void markBonds(atom *src, atom *neighbors, int distcount, int max);
void fixupLongBondChains(atom *src, atom *neighbors, int cutoff);
void updateHydrogenInfo(atom *allAtoms, atomBins *bins);
int atomsClose(atom *a, atom *b, float probeRad);
int inRange(point3d *p, point3d *q, float lim);
float gapSize(point3d *p, point3d *q, float qrad);
atom* findTouchingAtoms(atom *src, atomBins *bins, float probeRad);
void connectAll(FILE *outf, atom *allAtoms, atomBins *bins,
   float probeRad, char *id, int maxBnds, int noBonds);
int checkHB(atom *src, atom *targ, int maskHB);
void writeConnectivity(FILE *outf, atom *src, atom *neighList,
   char *id, int noBonds);
#endif
