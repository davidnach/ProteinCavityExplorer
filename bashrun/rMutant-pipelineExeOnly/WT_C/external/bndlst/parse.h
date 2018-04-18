/* name: parse.h                                         */
/* author: J. Michael Word     date written: 2/29/96     */
/* purpose: build a parse tree from an input string.     */
/*          based on "Compilers", Aho, Sethi, and Ullman */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef PARSE_H
#define PARSE_H 1

/* flags for atom properties */
#define         MC_PROP (1 <<  0)
#define         SC_PROP (1 <<  1)
#define      ALPHA_PROP (1 <<  2)
#define       BETA_PROP (1 <<  3)
#define          N_PROP (1 <<  4)
#define          C_PROP (1 <<  5)
#define          O_PROP (1 <<  6)
#define          S_PROP (1 <<  7)
#define          P_PROP (1 <<  8)
#define          H_PROP (1 <<  9)
#define      WATER_PROP (1 << 10)
#define        HET_PROP (1 << 11)
#define       PROT_PROP (1 << 12)
#define        DNA_PROP (1 << 13)
#define      METAL_PROP (1 << 14)
#define     METHYL_PROP (1 << 15)
#define      DONOR_PROP (1 << 16)
#define   ACCEPTOR_PROP (1 << 17)
#define   AROMATIC_PROP (1 << 18)
#define   CH_DONOR_PROP (1 << 19)
#define TEST_ACCEPT_ANGLE_PROP (1 << 20)

#define   NEGATIVE_PROP (1 << 21)
#define   POSITIVE_PROP (1 << 22)

#define   RHPHOBIC_PROP (1 << 23)
#define   RHPHILIC_PROP (1 << 24)
#define   RCHARGED_PROP (1 << 25)

#define   MAYBECHG_PROP (1 << 26)
#define CHECK_ENDS_PROP (1 << 27)

#endif
