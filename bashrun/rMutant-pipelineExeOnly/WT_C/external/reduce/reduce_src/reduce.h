// Name: reduce.C
// Author: J. Michael Word
// Date Written: 7/15/97
// Purpose: add hydrogens to a Protein Data Bank file and
//          writing a new Protein Data Bank file to standard output
//
// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999-2003 J. Michael Word
// **************************************************************
//
//  reduceChanges now contains the CHANGELOG or history info
//
 
#if defined(_MSC_VER)
#pragma warning(disable:4786) 
#pragma warning(disable:4305) 
#pragma warning(disable:4800) 
#endif

#include <iostream>
#include <fstream>

#ifdef OLD_STD_HDRS
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#else
#include <cstdlib>
#include <cstring>
#include <cctype>
#endif

#ifndef NOSYSSTATS
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif
#include <vector>
#include <list>
#include "utility.h"
#include "CTab.h"
#include "StdResH.h"
#include "ResBlk.h"
#include "pdb++.h"
#include "PDBrec.h"
#include "AtomPositions.h"

#define DIRECTORY_SEP_CHAR '/'

int ReturnCodeGlobal = 0;
#define ABANDONED_RC 1

using namespace std;

bool Verbose                  = FALSE;    // do we write processing notes to
					  // stdout?
bool KeepConnections          = TRUE;
bool StandardizeRHBondLengths = TRUE;
bool ProcessConnHydOnHets     = TRUE;
bool RemoveHydrogens          = FALSE;
bool BuildHisHydrogens        = FALSE;
bool SaveOHetcHydrogens       = TRUE;
bool UseXplorNames            = FALSE;
bool UseOldNames	      = FALSE; 
bool BackBoneModel	      = FALSE; 
bool DemandRotAllMethyls      = FALSE;
bool RotExistingOH            = FALSE;
bool DemandRotNH3             = TRUE;
bool DemandRotExisting        = FALSE;
bool DemandFlipAllHNQs        = FALSE;
bool DoOnlyAltA               = TRUE;
bool OKProcessMetMe           = TRUE;
bool OKtoAdjust               = TRUE;
bool ShowCliqueTicks          = TRUE;
bool ShowOrientScore          = FALSE;

int MinNTermResNo     = 1;    // how high can a resno be for n-term?
int ModelToProcess    = 1;    // which model to work on, >0 is a model
                              // to work on 041113
int ModelSpecified    = 0;    // commandline model specified 041113
int ModelNext         = 0;    // next model to process 041113
int ModelActive       = 0;    // found the next model and working on
			      // it 041113
int NBondCutoff       = 3;    // how many bonds away do we drop?
int ExhaustiveLimit   = 600;  // time limit, in seconds, to spend in
			      // brute force enumeration for a single
			      // clique
float ProbeRadius     = 0.0;  // how big is the probe in VDW
			      // calculations?
float VdwDotDensity   =16.0;  // how many dots per sq Angstroms in VDW
			      // calculations?
float OccupancyCutoff = 0.01; // lowest occupancy considered when
			      // determining score
float WaterBcutoff    =40.0;  // limit for water B values
float WaterOCCcutoff  = 0.66; // limit for water occupancy
float PenaltyMagnitude= 1.00; // score bias towards original
			      // orientation (changed from 0.0 in
			      // 2.13.0)
float MinRegHBgap     = 0.6;  // Hbonds with greater gaps start to
			      // bump
float MinChargedHBgap = 0.8;  // charged Hbonds start to bump at this
			      // point
float BadBumpGapCut   = 0.4;  // bump is bad if >= than this
float NonMetalBumpBias= 0.125;// bumps if H closer than atom radius,
			      // plus this
float MetalBumpBias   = 0.865;// ditto, for metals

string OFile; // if file exists, given orientations forced
bool UseSEGIDtoChainMap = FALSE; // if true, override some chain ids
char * inputPDBFile = NULL;

#ifndef HET_DICTIONARY
//#define HET_DICTIONARY "reduce_het_dict.txt"
#define HET_DICTIONARY "reduce_wwPDB_het_dict.txt"
#endif
string DBfilename( HET_DICTIONARY );

enum ConnType {NTERM_RES, CONNECTED_RES, FRAGMENT_RES};

struct SummaryStats {
  SummaryStats() : _H_found(0),        _H_HET_found(0),
		   _H_removed(0),      _H_HET_removed(0),
		   _H_added(0),        _H_HET_added(0),
		   _H_standardized(0), _H_HET_standardized(0),
		   _num_atoms(0),      _conect(0),
		   _num_adj(0),        _num_renamed(0) {}
  
  int _H_found, _H_HET_found, _H_removed, _H_HET_removed, _H_added;
  int _H_HET_added, _H_standardized, _H_HET_standardized, _num_atoms;
  int _conect, _num_adj, _num_renamed;
};

SummaryStats Tally;

char* parseCommandLine(int argc, char **argv);
void processPDBfile(istream& ifs, char *pdbFile, ostream& ofs);
void establishHetDictionaryFileName(void);
istream& inputRecords(istream& is, list<PDBrec*>& records);
ostream& outputRecords(ostream& os, const list<PDBrec*>& records);
void outputRecordsToFile(const list<PDBrec*>& records);
void dropHydrogens(list<PDBrec*>& records);
void reduceList(CTab& db, list<PDBrec*>& records,
		AtomPositions& xyz, vector<string>& fixNotes);
void scanAndGroupRecords(list<PDBrec*>& rlst, AtomPositions& xyz,
			 list<PDBrec*>::iterator& startAtoms);
void renumberAndReconnect(list<PDBrec*>& rlst);
void renumberAtoms(list<PDBrec*>& rlst);
void analyzeRes(CTab& db, ResBlk* pb, ResBlk* cb, ResBlk* nb,
		AtomPositions& xyz, list<PDBrec*>& waters, 
		vector<string>& fixNotes, list<PDBrec*>& rlst);
bool isConnected(ResBlk* a, ResBlk* b);
void genHydrogens(const atomPlacementPlan& pp, ResBlk& theRes, 
		  bool o2prime, AtomPositions& xyz, 
		  list<char>& resAlts, vector<string>& fixNotes, 
		  list<PDBrec*>& rlst);
void noteWaterInfo(ResBlk& r, list<PDBrec*>& waters);
void findAndStandardize(const char* name, const char* resname, 
			ResBlk& theRes);
void stdBondLen(float dist, PDBrec& ourHydrogen, 
		list<PDBrec*>& firstAtoms,
		const ElementInfo& e);
void fixupHeavyAtomElementType(ResBlk& theRes, CTab& hetDB);
void noteAddedInTally(const PDBrec& theH);
void noteStdInTally(const PDBrec& theH);
void noteRemovedInTally(const PDBrec& theH);
int findConnAtoms(const atomPlacementPlan& pp, ResBlk& theRes, 
		  char hac, PDBrec& r0atom, PDBrec& r1atom, 
		  PDBrec& r2atom);
bool okToPlaceHydHere(const PDBrec& theHatom,
		      const atomPlacementPlan& pp, 
		      const PDBrec& a1, const PDBrec& a2,
		      AtomPositions& xyz, bool& doNotAdjustSC,
		      vector<string>& fixNotes);
void recordSkipInfo(bool skipH, vector<string>& fixNotes,
		    const PDBrec& theHatom, const PDBrec& heavyAtom,
		    list<PDBrec*>& nearr, const char * msg);
void addHBondsToPDBFile(char *fileName);
