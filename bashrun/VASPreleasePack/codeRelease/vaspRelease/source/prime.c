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
 * File: prime.c
 *       prime code for set library
 *
 * Written by 
 *       Andrew M. Ladd <aladd@cs.rice.edu>
 * Modified by 
 *       Brian Y. Chen <chen@lehigh.edu>
 *
 * WWW URL: http://cse.lehigh.edu/~chen/
 * Email: chen@lehigh.edu
 * Documentation can be found here: 
 * http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
 *
 ***************************************************************************/

#include "defs.h"
#include "set.h"

#define MAX_NUMBER 65536
#define MAX_NUM_PRIMES 65536

static int erato_table[MAX_NUM_PRIMES];
static int total_primes = 0;

static void 
fill_table()
{
  int table[MAX_NUMBER];
  int i;
  
  for (i=0; i<MAX_NUMBER; i++) table[i] = 0;

  for (i=2; i<MAX_NUMBER; i++) 
    if (table[i] == 0) {
      int j;
      
      for (j=2*i; j<MAX_NUMBER; j+=i) table[j] = 1;
    }
  
  for (i=2; i<MAX_NUMBER; i++)
    if (table[i] == 0) erato_table[total_primes++] = i;
}


boolean_t_new 
is_prime(int n)
{
  int i;

  if (total_primes == 0) fill_table();

  for (i=0; i<total_primes; i++) {
    if (n % erato_table[i] == 0) return FALSE;
    if (erato_table[i] * erato_table[i] > n) return TRUE;
  }

  return TRUE;
}

int 
next_prime(int n)
{
	int trueStuff = 1;	///to elim stupid compiler warning
	int result = 0;

	while(trueStuff){
		if (is_prime(n)){
			result = n;
			break;
		}
		else n++;
	}

	return n;
}





/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/




