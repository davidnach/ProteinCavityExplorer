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
 * File: defs.h
 *       typedefinitions for set object
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

#ifndef ____RGL___DEFS___H___
#define ____RGL___DEFS___H___

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus 
#define BEGIN_C_DECL extern "C" {
#define END_C_DECL }
#else
#define BEGIN_C_DECL 
#define END_C_DECL 
#endif

BEGIN_C_DECL

/**
 * A more descriptive type for a general pointer.
 */

typedef void * ptr_t;

/**
 * A more descriptive type for a byte (8 bits).
 */

typedef unsigned char byte_t;

/**
 * A more descriptive type for a word (32 bits).
 */

typedef unsigned int word_t;

/**
 * A more descriptive type for a string.
 */

typedef char * string_t;

/**
 * A more descriptive type for a chunk of bytes.
 */

typedef byte_t * chunk_t;

/**
 * A more descriptive type for a boolean type.
 */

typedef int boolean_t_new;

#define TRUE 1
#define FALSE 0

#define NIL 0x0

END_C_DECL

#endif



/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/


