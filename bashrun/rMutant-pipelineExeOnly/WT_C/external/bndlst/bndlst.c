/* name: bndlst.c                        */
/* author: J. Michael Word               */
/* date written: 4/9/99                  */
/* purpose: output connectivity of PDB   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

/* modifications:                                                */
/* 4/09/99 - JMW - Modified from atvol V1.2.1                    */
/* 4/12/99 - JMW - Updated atom count for phantom water hydrogens*/
/*                 modified the output format, MSE not a HET,    */
/*                 select alt-conf code                          */
/* 5/31/00 - JMW - Added modified bases used in tRNAs            */
/*                 Included new properties TEST_ACCEPT_ANGLE and */
/*                 CH_DONOR                                      */
/*10/31/00 - JMW - Added O2* & O2' to NAbackboneList (oversite)  */
/*11/ 6/00 - JMW - Changed memory.h to stdlib.h                  */
/*11/19/00 - JMW - Updated select.c base list, added -nofacehb & */
/*                 included dna/rna (b) or protein (p) properties*/
/*11/27/00 - JMW - updated properties of odd bases in select.c   */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "bndlst.h"

static char *versionString =
	 "bndlst: version 1.6 11/27/00, Copyright 1999,2000, J. Michael Word";
static char *shortVersionStr = "Bndlst V1.6 (11/27/00)";

static int OnlyBonded= TRUE;  /* global controling output          */
static int ImplicitH = FALSE; /* global controling vdw radii       */
static int UseHParent= TRUE;  /* global controling H atom cov bond */
static int UseStdBond= TRUE;  /* global flag for Standard Bonding  */
static int HB2aromFace=TRUE;  /* global controling aromatic face Hbonding */
static int Verbose   = TRUE;  /* global flag for output messages   */
static int DebugLevel= 0;     /* global controling debugging information */

static float RadScaleFactor       =1.0; /* global VDW radius scale Factor */
static float CORadScale   =(1.65/1.75); /* global VDW radius scale Factor for C=O */

#define DEFAULTALT ('A')

static char TheAltCode = DEFAULTALT;/* which alt-conf code to select */
static long AtomCounter = 0; /* unique identifier for atoms */

#define USEPOLARH TRUE

#define QUARTER_PI (atan2(1.0,1.0))

#define MAX_BONDS (4)

#define ATOMS_TOO_CLOSE (0.000001)

#define TRACEATOM(where, txt, atptr, endtxt) \
(fprintf((where), "%s{%s%c %s%d%c}[%d]%s", (txt),\
		(atptr)->atomname, (atptr)->altConf,\
		(atptr)->resname, (atptr)->resid,\
		(atptr)->resInsCode, (atptr)->rescnt, (endtxt)))

int main(int argc, char **argv) {
   region boundingBox;
   atom *allAtoms   = NULL;
   atomBins *bins   = NULL;
   char *idString   = NULL;
   int  maxBonds    = 0;
   int   noBonds    = 0;
   int  modelNum    = 0;
   float probeRad   = 0.0;

   AtomCounter = 0;

   initStdConnTable();
   initalizeAtomTbl();

   allAtoms = processCommandline(argc, argv, &boundingBox,
      &idString, &maxBonds, &noBonds, &modelNum, &probeRad);

   if (allAtoms == NULL) {
      fprintf(stderr, "ERROR: No atom data in input.\n");
   }
   else {
      fprintf(stdout, "# output from %s\n", shortVersionStr);
      fprintf(stdout, "# max bonds = %d, cutoff = %d\n", maxBonds, noBonds);
      fprintf(stdout, "# probe radius = %g\n", probeRad);
      fprintf(stdout, "# radius scale factor = %g\n", RadScaleFactor);
      fprintf(stdout, "# C=O carbon VDW scaled by %.3f to a radius of %g A\n",
			CORadScale, CORadScale*1.75);
      fprintf(stdout, "# alt-conf selector = '%c'\n", TheAltCode);
      fprintf(stdout,"#type:id:  atomname  :atnum:    x   :    y   :    z   "
                     ":occ :bval:rad:ele:prop:n:dist:extras\n");

      if (Verbose) {
	 if (modelNum != 0) {
	    fprintf(stderr, "Model number %d selected\n", modelNum);
	 }
      }

      bins = binAtoms(allAtoms, &boundingBox, probeRad);

      if (! ImplicitH) {
	 updateHydrogenInfo(allAtoms, bins);
      }

      reverseAtomList(&allAtoms);

      connectAll(stdout, allAtoms, bins, probeRad, idString,
	 maxBonds, noBonds);
   }
   return 0;
}

void bndlsthelp() {
   fprintf(stderr, "\n");
   fprintf(stderr, "Syntax: bndlst input.pdb\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "Prints an entry for each atom followed by its neighbors\n");
   fprintf(stderr, "ending with an \"eol\" record. The last column on an\n");
   fprintf(stderr, "atom record contains the number of neighbors. The last\n");
   fprintf(stderr, "columns for a neighbor contain the bonding distance, \n");
   fprintf(stderr, "the Cartesian distance and notes on disulfides, etc..\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "Flags:\n");
   fprintf(stderr, "  -Name \"xxx\"  identifier for output table\n");
   fprintf(stderr, "  -Model#      model to sample (default 0)\n");
   fprintf(stderr, "  -ALTa        alternate conformation selected (default '%c')\n", DEFAULTALT);
   fprintf(stderr, "\n");
   fprintf(stderr, "  -NEARby      display all nearby atoms\n");
   fprintf(stderr, "  -Bonds#      how many bonds to trace (default %d)\n", MAX_BONDS);
   fprintf(stderr, "  -CUToff#     min num of bonds away (default 0)\n");
   fprintf(stderr, "  -Radius#.#   radius of the contact probe (default 0.0)\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "  -VDWscale#.# scale Van der Waals radii (default 1.0)\n");
   fprintf(stderr, "  -COSCale#.#  scale C=O carbon Van der Waals radii (default 0.94)\n");
   fprintf(stderr, "  -NOSTDbonds  do not consult the standard bonding table\n");
   fprintf(stderr, "  -NOHParent   do not use standard hydrogen-parent table\n");
   fprintf(stderr, "  -NOFACEhbond do not identify HBonds to aromatic faces\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "  -VErbose     verbose mode (default)\n");
   fprintf(stderr, "  -Quiet       quiet mode\n");
   fprintf(stderr, "  -Help        display this message\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "Atom properties listed in the output are coded as follows:\n");
   fprintf(stderr, "m mainchain, s sidechain, +/- charge, D donor, A acceptor\n");
   fprintf(stderr, "r aromatic ring, H het atom, w water, n nucleic acid,\n");
   fprintf(stderr, "p protein, M metalic.\n");
   fprintf(stderr, "\n%s\n", versionString);
   exit(0);
}

atom* processCommandline(int argc, char **argv, region *bbox,
   char **idString, int *maxBonds, int *noBonds,
   int *modelNum, float *probeRad) {
   int file = 0, i = 0, n = 0;
   char *p = NULL;
   FILE *inf = stdin;
   atom *atomlist = NULL;

   *idString   = "";
   *maxBonds   = MAX_BONDS;
   *noBonds    = 0;
   *modelNum   = 0;
   *probeRad   = 0.0;
   file = 0;

   for (i = 1; i < argc; i++) {
      p = argv[i];
      if (p[0] == '-') {
	 if ((file < 1) && (p[1] == '\0')) {
	    file = 1;
	    atomlist = loadAtoms(stdin, atomlist, bbox, *modelNum);
	 }
	 else if(compArgStr(p+1, "HELP", 1)){
	    bndlsthelp();
	 }
	 else if(compArgStr(p+1, "VERBOSE", 2)){
	    Verbose = TRUE;
	 }
	 else if(compArgStr(p+1, "QUIET", 1)){
	    Verbose = FALSE;
	 }
	 else if(compArgStr(p+1, "NEARby", 4)){
	    OnlyBonded = FALSE;
	 }
         else if(compArgStr(p+1, "NOFACEHBOND", 6)){
            HB2aromFace = FALSE;
         }
	 else if(n = compArgStr(p+1, "Radius", 1)){
	    *probeRad = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "VDWSCALE", 3)){
	    RadScaleFactor = parseReal(p, n+1, 10);
	    if (RadScaleFactor < 0.0001
	    ||  RadScaleFactor > 1000.0) {
	       fprintf(stderr, "ERROR: invalid scale factor -vdw%g\n",
			RadScaleFactor);
	       exit(0);
	    }
	 }
	 else if(compArgStr(p+1, "NAME", 1)){
	    if (++i < argc) {
	       *idString = argv[i];
	    }
	    else {
	       fprintf(stderr, "ERROR: no characters after -Name flag\n");
	       exit(0);
	    }
	 }
	 else if(n = compArgStr(p+1, "ALT", 3)){
	    TheAltCode = p[n+1];
	    if (TheAltCode == '\0') {
	       fprintf(stderr, "ERROR: truncated alt-conf selector, using blank\n");
	       TheAltCode = ' ';
	    }
	 }
	 else if(n = compArgStr(p+1, "COSCALE", 4)){
	    CORadScale = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "Bonds", 1)){
	    *maxBonds = parseInteger(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "CUToff", 3)){
	    *noBonds = parseInteger(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "Model", 1)){
	    *modelNum = parseInteger(p, n+1, 10);
	 }
	 else if(compArgStr(p+1, "STDBONDS", 7)){
	    UseStdBond = TRUE;
	 }
	 else if(compArgStr(p+1, "NOSTDbonds", 5)){
	    UseStdBond = FALSE;
	 }
	 else if(compArgStr(p+1, "HParent", 2)){
	    UseHParent = TRUE;
	 }
	 else if(compArgStr(p+1, "NOHParent", 4)){
	    UseHParent = FALSE;
	 }
	 else if(compArgStr(p+1, "IMPlicit", 3)){
	    ImplicitH = TRUE;
	 }
	 else if(n = compArgStr(p+1, "DEBUG", 5)){
	    DebugLevel = parseInteger(p, n+1, 10);
	 }
	 else {
	    fprintf(stderr, "ERROR: unrecognized flag, %s\n", p);
	 }
      }
      else if (file < 1) {
	 file = 1;
	 inf = fopen(p, "r");
	 if (inf) {
	    atomlist = loadAtoms(inf, atomlist, bbox, *modelNum);
	    fclose(inf);
	 }
	 else {
	    fprintf(stderr, "ERROR: could not open file %d: %s\n",
		     file, p);
	 }
      }
      else {
	 fprintf(stderr, "ERROR: duplicate filename, %s\n", p);
      }
   }
   if (atomlist == NULL) {
      bndlsthelp();
   }
   return atomlist;
}

void reverseAtomList(atom **allAtoms) {
   atom *a = NULL, *r = NULL, *nxt = NULL;

   a = *allAtoms;
   r = NULL;
   while (a) {
      nxt = a->next;
      a->next = r;
      r = a;
      a = nxt;
   }
   *allAtoms = r;
}

/* data structure used to determine if N and O are at chain ends */
typedef struct {
   atom *ambigN[8];    /* first residue's ambiguous N in this chain [0-3]*/
                       /* along with the last residue's Ns [4-7]         */
   atom *ambigO[8]; /*  last residue's ambiguous Os in this chain [0-7]  */
   int res_mc_oxy_cnt;       /* how many mainchain Oxygens on last res?  */
   int first, last;          /* where do the chain ids break?            */
   int Ntmarkers, Ctmarkers; /* where were N and C term markers observed? */
   int hisNHcount[4];
} chainEndData_t;

void initEndData(chainEndData_t *ed) {
   int i;
   for (i = 0; i < 8; i++) { ed->ambigN[i] = NULL; ed->ambigO[i] = NULL; }
   ed->res_mc_oxy_cnt = 0;
   ed->first = 0;
   ed->last  = 0;
   ed->Ntmarkers = 0;
   ed->Ctmarkers = 0;
   for (i = 0; i < 4; i++) { ed->hisNHcount[i] = 0; }
}

/* His ring atoms start out with a positive charge.                       */
/* Here we remove the charge for His residues without two ring NH protons.*/
void hisNHcheck(chainEndData_t *ed, atom *atomlist, int rescnt) {
   int i;
   atom *rat;

   if (ed->hisNHcount[0]) { /* the zeroth counter flags: Is this a HIS? */

      for(rat = atomlist; rat; rat = rat->next) {
	 if ((rat->rescnt == rescnt) && (rat->props & POSITIVE_PROP)) {
	    char altc = toupper(rat->altConf);
	    int ialt = 0;

	    if ((altc==' ')||(altc=='A')||(altc=='1')) { ialt = 1; }
	    else if ((altc=='B')||(altc=='2'))         { ialt = 2; }
	    else if ((altc=='C')||(altc=='3'))         { ialt = 3; }
	    
	    if (ialt > 0) {
	       if (ed->hisNHcount[ialt] < 2) { /* his not positive */
		  rat->props &= ~POSITIVE_PROP;
	       }
	       else {  /* his positive but not an acceptor */
		  rat->props &= ~ACCEPTOR_PROP;
	       }
	    }
	 }
      }

      for (i = 0; i < 4; i++) { ed->hisNHcount[i] = 0; }
   }
}

/* called after each residue is read */
void resCheck(chainEndData_t *ed, atom *atomlist, int rescnt) {
    hisNHcheck(ed, atomlist, rescnt);
}

/* called when a new residue begins and at the end of all processing */
/* to process the previous residue                                   */
void CtermCheck(chainEndData_t *ed, int rescnt, int isChainEnd) {
   int i;
   /* avoid processing the first time 'cause there ain't no chain pending */

   if (isChainEnd && rescnt) {
      ed->last = rescnt; /* last residue */

      /* see if we can put the pieces together to determine end charge */

      for (i = 0; i < 4; i++) { /* only array[0-3] contains first Ns */
	 if (ed->ambigN[i]) {
	    if (ed->ambigN[i]->rescnt == ed->first
	    &&  ed->Ntmarkers         == ed->first) {
	       ed->ambigN[i]->props |= POSITIVE_PROP;
	    }
	 }
	 else { break; }
      }
      for (i = 0; i < 8; i++) { /* array[0-7] contains last Os */
	 if (ed->ambigO[i]) {
	    if (ed->ambigO[i]->rescnt == ed->last
	    &&  ed->Ctmarkers         == ed->last) {
	       ed->ambigO[i]->props |= NEGATIVE_PROP;
	    }
	 }
	 else { break; }
      }

      initEndData(ed); /* reset the end data record */
   }
}

/* called when we have an indicator that residue is at N terminus */
/* to look for ambigN and mark charged */
void noticedNt(chainEndData_t *ed, int rescnt) {
   int i;
      for (i = 4; i < 8; i++) { /* only array[4-7] contains last Ns */
	 if (ed->ambigN[i]) {
	    if (ed->ambigN[i]->rescnt == rescnt) {
	       ed->ambigN[i]->props |= POSITIVE_PROP;
	    }
	 }
	 else { break; }
      }
}
/* called when we have an indicator that residue is at C terminus */
/* to look for ambigN and mark charged */
void noticedCt(chainEndData_t *ed, int rescnt) {
   int i;
      for (i = 0; i < 8; i++) { /* array[0-7] contains last Os */
	 if (ed->ambigO[i]) {
	    if (ed->ambigO[i]->rescnt == rescnt) {
	       ed->ambigO[i]->props |= NEGATIVE_PROP;
	    }
	 }
	 else { break; }
      }
}

/* called whenever a new residue begins */
void NtermCheck(chainEndData_t *ed, int rescnt, int isChainEnd) {
   int i;
   
   ed->res_mc_oxy_cnt = 0; /* reset Oxygen data for each residue */
   for (i = 0; i < 8; i++) { ed->ambigO[i] = NULL; }

   /* we have reserved the top 4 nitrogen slots for most recent res */
   /* so we have to clear this here                                 */
   for (i = 4; i < 8; i++) { ed->ambigN[i] = NULL; }
   
   if (isChainEnd) {
      ed->first = rescnt; /* first residue */
   }
}
/* called for each atom read in */
void ProcessResInfo(chainEndData_t *ed, atom *a) {
   int i;

   if (strstr("HIS:his", a->resname)) {
      ed->hisNHcount[0] = 1; /* indicates this is a HIS */
      if ( strstr(" HD1: hd1: HE2: he2", a->atomname)) {
	 char altc = toupper(a->altConf);
	 if ((altc==' ')||(altc=='A')||(altc=='1')) {
	    ed->hisNHcount[1]++;
	 }
	 else if ((altc=='B')||(altc=='2')) {
	    ed->hisNHcount[2]++;
	 }
	 else if ((altc=='C')||(altc=='3')) {
	    ed->hisNHcount[3]++;
	 }
      }
   }

   /* look for indicators that we are at the end of a chain */

   if (a->props & CHECK_ENDS_PROP) {
      if (isHatom(a->elem)) {
	 if (! ed->Ntmarkers) { ed->Ntmarkers = a->rescnt; }
	 noticedNt(ed, a->rescnt);
      }
      else if ( (a->elem == atomO) &&
	       ((a->altConf==' ')||(a->altConf=='A')||(a->altConf=='1')) ){
	 ed->res_mc_oxy_cnt++;
	 if (ed->res_mc_oxy_cnt == 2) {
	    ed->Ctmarkers = a->rescnt;
	    noticedCt(ed, a->rescnt);
	 }
      }
      else if (strstr(a->atomname, "OXT")) { /* check SUBSET of name */
	 ed->Ctmarkers = a->rescnt;
	 noticedCt(ed, a->rescnt);
      }
   }

   /* save pointers to the first ambigNs and the last ambigOs and Ns */
   /* we save multiple atoms to both find both possible mc oxygens   */
   /* and to handle multiple conformations (we can do 4 max)         */
   /* All atoms must have the same residue counter                   */

   if (a->props & MAYBECHG_PROP) {
      if (!strcmp(a->atomname, " N  ")) {
	 for (i = 0; i < 4; i++) {/* first Ns[0-3] (cleared at chain ends)*/
	    if (ed->ambigN[i] == NULL) {
	       if ((i == 0) || (ed->ambigN[i-1]->rescnt == a->rescnt)) {
		  ed->ambigN[i] = a;
	       }
	       break;
	    }
	 }
	 for (i = 4; i < 8; i++) {   /* last Ns[4-7] (cleared each res)*/
	    if (ed->ambigN[i] == NULL) {
	       if ((i == 4) || (ed->ambigN[i-1]->rescnt == a->rescnt)) {
		  ed->ambigN[i] = a;
	       }
	       break;
	    }
	 }
      }
      if (!strcmp(a->atomname, " O  ")) { /* last Os[0-7]  (cleared each res)*/
	 for (i = 0; i < 8; i++) {
	    if (ed->ambigO[i] == NULL) {
	       if ((i == 0)|| (ed->ambigN[i-1]->rescnt == a->rescnt)) {
		  ed->ambigO[i] = a;
	       }
	       break;
	    }
	 }
      }
   }

}
atom* loadAtoms(FILE *fp, atom *atomlist, region *boundingBox, int modelNum) {
   int previd = 0, rescnt = 0, model = 0;
   char *rec, prevInsCode = ' ', prevChain = '?';
   atom *a = NULL;
   chainEndData_t endData;

   initEndData(&endData);

   while(rec = getPDBrecord(fp)) {
      if (isTer(rec)) { prevChain = '?'; }
      if ((isAtom(rec) || isHet(rec)) && ! isPseudoAtom(rec) && (model == modelNum)) {
	 a = newAtom(rec, 1, model);
	 if (a) {
	    if (a->resid != previd || a->resInsCode != prevInsCode
	                           || a->chain      != prevChain) {

	       if (rescnt){
		  resCheck(&endData, atomlist, rescnt);
	       }
	       CtermCheck(&endData, rescnt, (a->chain != prevChain));

	       ++rescnt;

	       NtermCheck(&endData, rescnt, (a->chain != prevChain));
	       previd = a->resid;
	       prevInsCode = a->resInsCode;
	       prevChain = a->chain;
	    }
	    a->rescnt = rescnt;

	    ProcessResInfo(&endData, a);

	    if (ImplicitH && isHatom(a->elem)) {
	       free(a); /* filter out implicit hydrogens */
	       continue;
	    }

	    if (!atomlist) {
	       boundingBox->min = a->loc;
	       boundingBox->max = a->loc;
	    }
	    a->next = atomlist;
	    atomlist = a;

	    updateBoundingBox(&(a->loc), boundingBox);
	 }
      }
      if (isModel(rec)) {
	 model = parseModel(rec);
      }
   }
   if (rescnt) { resCheck(&endData, atomlist, rescnt); }
   CtermCheck(&endData, rescnt, TRUE);

   return atomlist;
}

atomBins* binAtoms(atom *allAtoms, region *boundingBox, float probeRad) {
   atom *a;
   atomBins *bins;

   bins = initBins(boundingBox,
	    2.0*(getMaxRadius(ImplicitH)+probeRad) + 0.2);
   if (bins) {
      for(a = allAtoms; a; a = a->next) {
	 addNeighbor(a, bins);
      }
   }
   /* note: addNeighbor() also called in updateHydrogenInfo() */
   return bins;
}

float getRadius(int at) {
   float rad = 0.0;

   if (ImplicitH) { rad = getImplRad(at); }
   else           { rad = getExplRad(at); }

   return rad*RadScaleFactor;
}

atom * newAtom(char *rec, int file, int model) {
   atom *a;

   a = (atom *)malloc(sizeof(atom));
   if (a) {
      a->next     = NULL;
      a->nextInBin= NULL;
      a->scratch  = NULL;
      a->mark = 0;    /* used to mark bonds */
      a->atomCount = ++AtomCounter;
      a->flags = 0;
      a->props = 0;
      parseResidueName(rec, a->resname);
      parseAtomName(rec, a->atomname);
      a->file = file;
      a->model = model;
      a->chain = parseChain(rec);
      a->resid = parseResidueNumber(rec);
      a->resInsCode = parseResidueInsertionCode(rec);
      a->rescnt = 0;
      a->altConf = parseAltLocCode(rec);
      if (strstr(":ASX:GLX:ASN:GLN:", a->resname) /* special case treats undecided */
      && (a->atomname[1] == 'A')) {             /* as an oxygen */
	 a->elem = identifyAtom(" O  ");
	 fprintf(stderr, "WARNING: atom %s will be treated as oxygen\n", a->atomname);
      }
      else { /* normal case */
	 a->elem = identifyAtom(a->atomname);
      }

      if (ImplicitH && isHatom(a->elem)) { free(a); return NULL; }

      a->bondedto = NULL;
      if (isHatom(a->elem)) {
	 a->atomclass = -1; /* to be specified later */
	 if (UseHParent) {
	    a->bondedto = setHydrogenParentName(a->resname, a->atomname);
	 }
      }
      else {
	 a->atomclass = a->elem;
	 if (UseStdBond) {
	    a->bondedto = setMainchainBonding(a->resname, a->atomname);
	 }
      }

      setProperties(a, isHet(rec), HB2aromFace);

      parseXYZ(rec, &(a->loc));

      a->radius = getRadius(a->elem);
      a->covRad = getCovRad(a->elem);

   /* note: H DONOR_PROP assignment done later in updateHydrogenInfo() */

      a->occ  = parseOccupancy(rec);
      a->bval = parseTempFactor(rec);
      parseSegID(rec, a->segid);

      if (isCarbonylAtom(a)) {
	 a->radius = CORadScale*getRadius(a->elem);
      }
   }
   else {
      fprintf(stderr, "ERROR: could not allocate space for new atom\n");
   }
   return a;
}

void markBonds(atom *src, atom *neighbors, int distcount, int max) {
   atom *targ, *theHatom, *theOtherAtom;
   int nearpt, tooclose, isAanH, isBanH, letItBond = FALSE;

   for(targ = neighbors; targ; targ = targ->scratch) {
      if (targ->mark == 0 || targ->mark > distcount) {

	 /* note: when we are relying on the StdBond table   */
	 /*       to figure out which bonds are allowed      */
	 /*       we can be more generous with the distance  */
	 /*       cuttofs to permit more oddball distortions */

	 nearpt = inRange(&(src->loc), &(targ->loc),
		     src->covRad + targ->covRad + COVRADFUDGE
		     + (UseStdBond ? COVRADFUDGE : 0.0) );

	 tooclose = inRange(&(src->loc), &(targ->loc),
			   src->covRad + targ->covRad - 0.6
			   - (UseStdBond ? COVRADFUDGE : 0.0) );

	 /* conditions for allowing a bond to be formed... */

	 if (nearpt && (! tooclose) ) {
	    isAanH = isHatom(src->elem);
	    isBanH = isHatom(targ->elem);

#ifdef OLD_BONDING_CALC
	    if ( (! isAanH && ! isBanH)             /* neither atom is a hydrogen  */
	    ||  ((! isAanH || ! isBanH)             /* or one is an H and both are */
	      && ATOMS_IN_THE_SAME_RES(src, targ)   /*    in the same res with parent name*/
	      && ((theHatom->bondedto == NULL)      /*     either unknown or as expected  */
	         || strstr(theHatom->bondedto, theOtherAtom->atomname) )) ) {
	       targ->mark = distcount;
	       if (distcount < max) {
		  markBonds(targ, neighbors, distcount+1, max);
	       }
	    }
#else
	    letItBond = FALSE; /* we start out skeptical */

	    if ((! isAanH) && (! isBanH)) { /* neither atom is a hydrogen - both are heavy atoms */

	       if (UseStdBond) { /* if we are allowing only "standard" bonding */

		  if (! ATOMS_IN_THE_SAME_RES(src, targ)) { /* and are in diff res */

		     if ((src->props & MC_PROP) && (targ->props & MC_PROP)) {

			/* both are mc in diff res - only specific bonding patterns allowed */

			letItBond =
			       (!strcmp(" N  ", src->atomname) && !strcmp(" C  ",targ->atomname))
			   ||  (!strcmp(" N  ",targ->atomname) && !strcmp(" C  ", src->atomname))
			   ||  (!strcmp(" P  ", src->atomname) &&  strstr(" O3*: O3'", targ->atomname))
			   ||  (!strcmp(" P  ",targ->atomname) &&  strstr(" O3*: O3'",  src->atomname));
		     }
		     else if (!(src->props & MC_PROP) && !(targ->props & MC_PROP)) {

			/* both are sc in diff res they must be CYS-CYS SS */

			letItBond = (!strcmp(src->resname,  "CYS") && !strcmp(src->atomname,  " SG ")
				 &&  !strcmp(targ->resname, "CYS") && !strcmp(targ->atomname, " SG "));
		     }
		     else { letItBond = FALSE; } /* sc/mc in diff res - bonding not possible */
		  }
		  else { /* both heavy atoms in same res */
#ifdef STRICT_CONN_TABLE
		     /* strict - heavy atom must be as exactally expected */

		     letItBond = ( (src->bondedto != NULL) && strstr(src->bondedto, targ->atomname) );
#else
		     /* heavy atom either unknown or as expected */

		     letItBond = ( (src->bondedto == NULL) || strstr(src->bondedto, targ->atomname) );
#endif
		  }
	       }
	       else { /* not using std bond table - so we let heavy atoms bond when they are close enough */
		  letItBond = TRUE;
	       }
	    }
	    else if ((! isAanH || ! isBanH)          /* only one atom is a hydrogen        */
	      && ATOMS_IN_THE_SAME_RES(src, targ)) { /* and both atoms are in the same res */

	       if (isAanH) { theHatom = src;  theOtherAtom = targ; }
	       else        { theHatom = targ; theOtherAtom = src;  }

	       if ( (theHatom->bondedto == NULL) /* heavy atom either unknown or as expected */
	       ||    strstr(theHatom->bondedto, theOtherAtom->atomname)
	       ||   (UseHParent == FALSE) ) {
		  letItBond = TRUE;
	       }
	    }
	    else { /* either both are hydrogens in the same residue - no bond possible */
	           /*        or one is a H and they are in diff res - no bond possible */
	       letItBond = FALSE;
	    }

	    if (letItBond) {

	       targ->mark = distcount;
	       if (distcount < max) {
		  markBonds(targ, neighbors, distcount+1, max);
	       }
#ifndef DUMP_BONDS
	       if (DebugLevel > 7) {
		  fprintf(stdout, "{%4.4s%c%3.3s %c%4d%c}P %8.3f%8.3f%8.3f\n",
			   src->atomname, src->altConf,
			   src->resname, src->chain,
			   src->resid, src->resInsCode,
			   src->loc.x, src->loc.y, src->loc.z);
		  fprintf(stdout, "{%4.4s%c%3.3s %c%4d%c}L %8.3f%8.3f%8.3f\n",
			   targ->atomname, targ->altConf,
			   targ->resname, targ->chain,
			   targ->resid, targ->resInsCode,
			   targ->loc.x, targ->loc.y, targ->loc.z);
	       }
#endif
	    }
#endif
	 }
      }
   }
}

/* fixupLongBondChains() - unmark and remove from the bonded set    */
/*                         atoms which are more than cutoff bonds   */
/*                         from the source when neither src or targ */
/*                         is a hydrogen                            */
void fixupLongBondChains(atom *src, atom *neighbors, int cutoff) {
   atom *targ;
   /* int nearpt; */

   if (! isHatom(src->elem)) {
      for(targ = neighbors; targ; targ = targ->scratch) {
	 if (targ->mark > cutoff && ! isHatom(targ->elem)) {
	    targ->mark = 0;
	 }
      }
   }
}

int dotType(atom *src, atom *atomList) {
   atom *a;
   int rslt = src->atomclass;

   if (src->atomclass < 0) {
      rslt = src->elem;

      if (isHatom(src->elem)) {
	 for(a = atomList; a; a = a->scratch) {
	    if (a->mark == 1 && ! isHatom(a->elem)
	                     && src->rescnt == a->rescnt) {
	       rslt = a->elem; break;
	    }
	 }
      }
      src->atomclass = rslt;
   }

   return rslt;
}

/* Make polor hydrogens HB donors.                          */
/* Waters without hydrogens are protonated with "phantoms". */

void updateHydrogenInfo(atom *allAtoms,atomBins *bins) {
   atom *src, *orig, *atomList, *a, *newH;
   int type, newHcount, i, nProx;

/* Information specific to making sure only one H2O Hydrogen is     */
/* generated which points into ring atoms on a given aromatic ring. */
   struct {atom  *a; float gap;} proximity[20];

   for(src = allAtoms; src; src = src->next) {
/* Note: we have to examine ALL the atoms at this point.        */

      if (isHatom(src->elem)) { /* fixup Hs in the input */

	 atomList = findTouchingAtoms(src, bins, 0.0);
	 if (atomList) {
	    markBonds(src, atomList, 1, 1);
	    type = dotType(src, atomList);

	    if ( type == atomN || type == atomO || type == atomS ) {

	       src->props |= DONOR_PROP; /* polar hydrogens are donors */

	       if (USEPOLARH) {
		  src->elem = atomHpolar; /* very important */
		  src->radius = getRadius(src->elem);
	       }
/* there are a more restrictive set of cutoffs in reduce: 0.7 & 40 */
#define H2OoccCUTTOFF   0.25
#define H2ObvalCUTTOFF 80.0
	       if ((FABS(src->occ) < H2OoccCUTTOFF
	         || src->bval >= H2ObvalCUTTOFF)
	       &&  (src->props & WATER_PROP)) {
		  src->elem = ignoreAtom;
		  src->radius = 0.0;
	       } /* ignore low occupancy water hydrogens */

	       for(a = atomList; a; a = a->scratch) {
	       /* since we found a hydrogen, then we know */
	       /* that the heavy atom is not the donor    */
	       /* so we can disable that property         */
		  if ((a->mark == 1) && ! isHatom(a->elem)
		  &&  (FABS(src->occ) >= H2OoccCUTTOFF)
		  &&  (src->bval < H2ObvalCUTTOFF)) {
		     a->props &= ~DONOR_PROP;
		  }
	       }
	    } /* end if (atom types) */
	 }
      }
	 /* fixup water oxygen */
      else if ((src->props & WATER_PROP) && (src->elem == atomO)) {

	 /*if connected H already found, the DONOR prop will not be set*/
	 /*if it is, we look further ahead in this residue for Hs      */
	 if (src->props & DONOR_PROP) {
	    for(a = src->next;
	        a && (src->rescnt == a->rescnt); a = a->next) {
	       if (  isHatom(a->elem)
	         && (src->altConf == a->altConf)
	         && (a->occ > 0.1)) {
		  src->props &= ~DONOR_PROP; /* turn off donor property */
	       }
	    }
	 }
	                                /* if donor flag still set...  */
	 if (src->props & DONOR_PROP) { /* water O without occupied Hs */
	    src->props |=  ACCEPTOR_PROP;
	    src->props &= ~DONOR_PROP;    /* acceptor, not donor */
	    orig = src; /* keep track of the original Oxygen atom */
	    atomList = findTouchingAtoms(src, bins, 0.25);
   /* the 0.25 comes from:
      OxygenVDW(1.4) + 2*probeRad(0.3) == HpolVDW(1.0) + OHBondLen(1.0)
      and we want a minimum overlap of 0.1 before we consider possible Hs
   */

	    /* for each nearby HB acceptor, keep only 1 per aromatic ring */
	    nProx = 0;
	    for(a = atomList; a; a = a->scratch) {
	       if (a->props & ACCEPTOR_PROP) {
		  float hbgap = gapSize(&(a->loc), &(orig->loc),
		                          a->radius + 2.0);

		  int foundAromRingAtom = FALSE;

		  if (a->props & AROMATIC_PROP) {
		     for (i=0; i < nProx; i++) {
			if (a->rescnt == proximity[i].a->rescnt
			     && proximity[i].a->props & AROMATIC_PROP) {

			   if (hbgap < proximity[i].gap) {
			      proximity[i].a = a;
			      proximity[i].gap = hbgap;
			   }
			   foundAromRingAtom = TRUE;
			   break;
			}
		     }
		  }
		  if ( (!foundAromRingAtom)
		  && nProx < sizeof(proximity)/sizeof(proximity[0])) {
		     proximity[nProx].a = a;
		     proximity[nProx].gap  = hbgap;		   
		     ++nProx;
		  }
	       } /* end if acceptor_prop */
	    }

	    /* build new HOH hydrogens in the direction of each acceptor */
	    newHcount = 0;
	    for (i=0; i < nProx; i++) {

	       a = proximity[i].a;

	       if (a->props & ACCEPTOR_PROP) {

		  /* Make a phantom H in the direction of each Acceptor */
		  newH = (atom *)malloc(sizeof(atom));
		  if (newH) {
		     point3d o2aVec;
		     double waterOHbondLen;

		     *newH = *orig;
		     newH->next = src->next; /* hook into list of atoms */
		     src->next  = newH;

		     v3sub(&(a->loc), &(orig->loc), &o2aVec);
#define BEST_HBOND_OVERLAP 0.6
		     waterOHbondLen = 1.0 +
			   MAX(-1.0, MIN(0.0, proximity[i].gap + BEST_HBOND_OVERLAP));
		     v3scale(&o2aVec, waterOHbondLen);
		     v3add(&(orig->loc), &o2aVec, &(newH->loc));

		     newH->nextInBin= NULL; /* added to bins below */
		     newH->scratch  = NULL;
		     newH->atomCount = ++AtomCounter;
		     newH->elem     = atomHOd;
		     newH->radius   = getRadius(newH->elem);
		     newH->covRad   = getCovRad(newH->elem);
		     strcpy(newH->atomname, " H? ");

		     newH->props |=  DONOR_PROP;
		     newH->props &= ~ACCEPTOR_PROP; /* donor, not acceptor */

		     addNeighbor(newH, bins); /* adding an atom not in input! */

if (DebugLevel>1) {
fprintf(stderr,
"HETATM%5d  H%-2d%c%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      new  H\n",
			   0, ++newHcount, newH->altConf,
			   newH->resname, newH->chain,
			   newH->resid, newH->resInsCode,
			   newH->loc.x, newH->loc.y, newH->loc.z,
			   newH->occ, newH->bval);
}

		     src = src->next; /* bypass new atom when we iterate */
		  }
	       } /* end if acceptor_prop */
	    } /* end for */

	 }
      } /* else if water */

      /* Remove any DONOR property for non-water polar atoms     */
      /* since the Hydrogen atom itself will be the DONOR.       */
      /* This was assigned in select.c to facilitate selections. */
      else if (!(src->props & WATER_PROP)
	    &&  (src->elem == atomN
	      || src->elem == atomO
	      || src->elem == atomS)) { src->props &= ~DONOR_PROP; }

   } /* end loop over all atoms */
}

int atomsClose(atom *a, atom *b, float probeRad) {
   int nearpt = FALSE;
   float lim, dsq;

   lim = a->radius + b->radius + probeRad + probeRad;

   dsq = v3distanceSq(&(a->loc), &(b->loc));

   /* if too close they must be the same atom... */

   if ((dsq > 0.001) && (dsq <= (lim*lim))) {
      nearpt = TRUE;
   }
   return nearpt;
}

int inRange(point3d *p, point3d *q, float lim) {
   return v3distanceSq(p, q) <= (lim*lim);
}

float gapSize(point3d *p, point3d *q, float qrad) {
   return v3distance(p, q) - qrad;
}

atom* findTouchingAtoms(atom *src, atomBins *bins, float probeRad){
   atom *a, *head = NULL;
   int i, j, k, nearpt, targcount = 0;

   if (src->ix < 1 || src->iy < 1 || src->iz < 1) { return NULL; }

   for(i = src->ix-1; i <= src->ix+1; i++) {
      for(j = src->iy-1; j <= src->iy+1; j++) {
	 for(k = src->iz-1; k <= src->iz+1; k++) {

	    for(a = bins->list[i][j][k]; a; a = a->nextInBin) {

	       if (( src->rescnt  == a->rescnt)  && /* alternate      */
	           ( src->altConf != a->altConf) && /* conformations  */
	           ( src->altConf != ' ')        && /* don't interact */
	           (   a->altConf != ' ')     ) { continue; }

	       nearpt = atomsClose(src, a, probeRad);
	       if (nearpt) {
		  a->mark = 0;       /* clear the bonding marker.   */
		  a->scratch = head; /* link using scratch pointers */
		  head = a;          /* with <a> at head of list    */

		  targcount++;
	       }
	    }
	 }
      }
   }
   src->mark = 0;       /* clear the src atom bonding marker.   */

   return (targcount < 1) ? NULL : head;
}

void connectAll(FILE *outf, atom *allAtoms, atomBins *bins,
	    float probeRad, char *id, int maxBonds, int noBonds) {
   atom *src, *neighList;

   for(src = allAtoms; src; src = src->next) {
      neighList = findTouchingAtoms(src, bins, probeRad);
      markBonds(src, neighList, 1, maxBonds);
      writeConnectivity(outf, src, neighList, id, noBonds);
   }
}

int checkHB(atom *src, atom *targ, int maskHB) {
   int bothCharged = FALSE, chargeComplement = FALSE;

   if ((src == NULL) || (targ == NULL)) { return 0; }

   if (( src->props & (NEGATIVE_PROP|POSITIVE_PROP))
   &&  (targ->props & (NEGATIVE_PROP|POSITIVE_PROP))) {
      bothCharged = TRUE;
      chargeComplement =
		     ((( src->props & POSITIVE_PROP) &&
		       (targ->props & NEGATIVE_PROP))
		  ||  (( src->props & NEGATIVE_PROP) &&
		       (targ->props & POSITIVE_PROP)));
   }

   return ((targ->props & maskHB)
	 && ((! bothCharged) || chargeComplement));
}

void writeConnectivity(FILE *outf, atom *src, atom *neighList,
                  char *id, int noBonds) {
   atom *a = NULL;
   point3d p;
   int i = 0, nneigh = 0, maskHB = 0, isaHBond = FALSE;
   /* int bothCharged = FALSE, chargeComplement = FALSE; */
   int countIt = FALSE;

   if (src->elem == ignoreAtom) { return; }

   if ( (TheAltCode   != ' ')
     && (src->altConf != ' ')
     && (src->altConf != TheAltCode)) { return; }

   maskHB = 0;
   if(src->props & DONOR_PROP) {
      maskHB = ACCEPTOR_PROP;
      if(src->props & ACCEPTOR_PROP) { /* ambig */
	 maskHB |= DONOR_PROP;
      }
   }
   else if (src->props & ACCEPTOR_PROP) {
      maskHB = DONOR_PROP;
   }

   fprintf(outf, "atom:%s", (id ? id : ""));
   fprintf(outf, ":%c%4d%c%-3.3s %4.4s%c:%5ld",
	 src->chain,src->resid, src->resInsCode, src->resname,
	 src->atomname, src->altConf, src->atomCount);

   fprintf(outf, ":%8.3f:%8.3f:%8.3f",src->loc.x,src->loc.y,src->loc.z);
   fprintf(outf, ":%4.2f:%4.1f",src->occ,src->bval);

   fprintf(outf, ":%4.2f:%s:", src->radius, getAtomName(src->elem));
   if (src->props & MC_PROP)       { fprintf(outf, "m"); }
   if (src->props & SC_PROP)       { fprintf(outf, "s"); }
   if (src->props & WATER_PROP)    { fprintf(outf, "w"); }
   if (src->props & POSITIVE_PROP) { fprintf(outf, "+"); }
   if (src->props & NEGATIVE_PROP) { fprintf(outf, "-"); }
   if (src->props & DONOR_PROP)    { fprintf(outf, "D"); }
   if (src->props & ACCEPTOR_PROP) { fprintf(outf, "A"); }
   if (src->props & AROMATIC_PROP) { fprintf(outf, "r"); }
   if (src->props & HET_PROP)      { fprintf(outf, "H"); }
   if (src->props & DNA_PROP)      { fprintf(outf, "n"); }
   if (src->props & PROT_PROP)     { fprintf(outf, "p"); }
   if (src->props & METAL_PROP)    { fprintf(outf, "M"); }

   nneigh = 0;
   for(a = neighList; a; a = a->scratch) {
      isaHBond = checkHB(src, a, maskHB);

      countIt = FALSE;
      if( (TheAltCode == ' ') ?
          ((src->altConf == ' ') || (a->altConf == ' ') || (src->altConf == a->altConf))
         :
          (((src->altConf == ' ') || (src->altConf == TheAltCode))
          && ((a->altConf == ' ') ||   (a->altConf == TheAltCode)))    )
      {

	 if (OnlyBonded) {
	    if ((a->mark == 1) || (isaHBond)) { countIt = TRUE; }
	 }
	 else {
	    if ((a->mark == 0) || (a->mark >= noBonds)) {
	       countIt = TRUE;
	    }
	 }
      }

      if (countIt) { nneigh++; }
   }
   fprintf(outf, ":%d:\n", nneigh);

   i = 0;
   for(a = neighList; a; a = a->scratch) {
      isaHBond = checkHB(src, a, maskHB);

      v3sub(&(a->loc), &(src->loc), &p);

      countIt = FALSE;
      if( (TheAltCode == ' ') ?
          ((src->altConf == ' ') || (a->altConf == ' ') || (src->altConf == a->altConf))
         :
          (((src->altConf == ' ') || (src->altConf == TheAltCode))
          && ((a->altConf == ' ') ||   (a->altConf == TheAltCode)))    )
      {

	 if (OnlyBonded) {
	    if ((a->mark == 1) || (isaHBond)) { countIt = TRUE; }
	 }
	 else {
	    if ((a->mark == 0) || (a->mark >= noBonds)) {
	       countIt = TRUE;
	    }
	 }
      }

      if (countIt) {
	 fprintf(outf, "%3d%c", ++i, ((a->mark == 1)?'b':(isaHBond?'h':' ')));
	 fprintf(outf, ":%s", (id ? id : ""));

	 fprintf(outf, ":%c%4d%c%-3.3s %4.4s%c:%5ld",
                           a->chain,a->resid, a->resInsCode, a->resname,
			   a->atomname, a->altConf, a->atomCount);
	 fprintf(outf, ":%8.3f:%8.3f:%8.3f",a->loc.x,a->loc.y,a->loc.z);
	 fprintf(outf, ":%4.2f:%4.1f",a->occ,a->bval);
	 fprintf(outf, ":%4.2f:%s:", a->radius, getAtomName(a->elem));

	 if (a->props & MC_PROP)       { fprintf(outf, "m"); }
	 if (a->props & SC_PROP)       { fprintf(outf, "s"); }
	 if (a->props & WATER_PROP)    { fprintf(outf, "w"); }
	 if (a->props & POSITIVE_PROP) { fprintf(outf, "+"); }
	 if (a->props & NEGATIVE_PROP) { fprintf(outf, "-"); }
	 if (a->props & DONOR_PROP)    { fprintf(outf, "D"); }
	 if (a->props & ACCEPTOR_PROP) { fprintf(outf, "A"); }
	 if (a->props & AROMATIC_PROP) { fprintf(outf, "r"); }
	 if (a->props & HET_PROP)      { fprintf(outf, "H"); }
	 if (a->props & DNA_PROP)      { fprintf(outf, "n"); }
	 if (a->props & PROT_PROP)     { fprintf(outf, "p"); }
	 if (a->props & METAL_PROP)    { fprintf(outf, "M"); }

	 if (a->mark > 0) {
	    fprintf(outf, ":%1d", a->mark);
	 }
	 else {
	    fprintf(outf, ":>");
	 }
	 fprintf(outf, ":%.3f:", v3length(&p));
	 
	 if ((src->elem == atomS) && (a->elem == atomS)) {
	    if (a->mark == 1) { fprintf(outf, " disulfide"); }
	 }
	 else if(((src->props & METAL_PROP) && !(a->props & METAL_PROP))
	     || (!(src->props & METAL_PROP) &&  (a->props & METAL_PROP))) {
	    if (a->mark == 1) { fprintf(outf, " metal-ligand"); }
	 }
	 else if(((src->props & HET_PROP) && !(a->props & HET_PROP))
	     || (!(src->props & HET_PROP) &&  (a->props & HET_PROP))) {
	    if (a->mark == 1) { fprintf(outf, " het-modified"); }
	 }
	 else if(((src->props & POSITIVE_PROP) && (a->props & NEGATIVE_PROP))
	     || ( (src->props & NEGATIVE_PROP) && (a->props & POSITIVE_PROP))) {
	    if (a->mark != 1) { fprintf(outf, " charge-complement"); }
	 }
	 else if(((src->props & POSITIVE_PROP) && (a->props & POSITIVE_PROP))
	     || ( (src->props & NEGATIVE_PROP) && (a->props & NEGATIVE_PROP))) {
	    if (a->mark != 1) { fprintf(outf, " charge-repulsion"); }
	 }

	 fprintf(outf, "\n");

      }
   }
   fprintf(outf, "eol :%s:\n", (id ? id : ""));
}
