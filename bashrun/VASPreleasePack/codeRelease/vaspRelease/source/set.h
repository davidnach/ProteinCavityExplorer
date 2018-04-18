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
 * File: set.h
 *       set library for fast lookup
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

#ifndef _____RGL_SET_H_____
#define _____RGL_SET_H_____

#include "defs.h"

typedef struct 
{
  int capacity;
  int size;
  
  int table_size;

  word_t properties;

  ptr_t header;
} set_head_t;

#define get_set_head(s) (set_head_t *) (((byte_t *) s) - sizeof(set_head_t))


/**
 * A set type which can also be used as an array.  The array should not 
 * written into.
 */

typedef int * set_t;

/**
 * An enumeration of set properties.
 */

typedef enum
  {
    SP_MAP = 0x01            /* Adds a map to the set. */
  } set_properties_t;

#define NO_SUCH_ELEMENT -1

BEGIN_C_DECL

/**
 * Allocates a set object.
 *
 * @param properties a union of set properties masks.
 */

set_t alloc_set(word_t properties);

/**
 * Deallocates a set object.
 */

void free_set(set_t s);

/**
 * Returns the number of elements placed in the set.
 */

int size_set(set_t s);

/**
 * Gets the set header (a pointer that can be freely assigned).
 */

ptr_t get_header_set(set_t s);

/**
 * Sets the set header (a pointer that can be freely assigned).
 */

void set_header_set(set_t s, ptr_t header);

/**
 * Puts an element into the set.
 * @param s a set.
 * @param id the 'id' of the element to put.
 * @return the new set - might be different.
 */

set_t put_set(set_t s, int id);

/**
 * Removes an element from the set.
 * @param s a set.
 * @param id the 'id' of the element to remove.
 * @return TRUE if an element was removed or FALSE if nothing removed.
 */

boolean_t_new remove_set(set_t s, int id);

/**
 * Swap element.
 * @param s a set.
 * @param i the first index.
 * @param j the second index.
 */

void swap_set(set_t s, int i, int j);

/**
 * Checks to see if a set contains a given element.
 * @param s a set.
 * @param id the 'id' of the element to be contained.
 * @return TRUE if the element in the set or FALSE otherwise.
 */

boolean_t_new contains_set(set_t s, int id);

/**
 * Returns the index of a given element in the set.
 * @param s a set.
 * @param id the 'id' of the element of a given element.
 * @return the index of the element or NO_SUCH_ELEMENT if no such element exists.
 */

int index_of_set(set_t s, int id);

/**
 * Associates an object with an element.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @param id the element to associate with.
 * @param obj the object to associate.
 * @return the new set - may be different.
 */

set_t associate_set(set_t s, int id, ptr_t obj);

/**
 * Returns the object mapped to a given element in the set.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @param id the element to check.
 * @return the object associated with element 'id' or NIL if not possible.
 */
 

ptr_t mapsto_set(set_t s, int id);

/**
 * The object stored at the ith index.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @param index the index of the element to use.
 * @return the object associated with the indexth element.
 */

ptr_t ith_map_set(set_t s, int index);

//ptr_t ith_map_set(set_t s, int index);

/**
 * Returns an array of objects associated in the set.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @return an array of objects.
 */

ptr_t * maptable_set(set_t s);

/**
 * Writes a set in human readable.
 * @param fp a file pointer.
 * @param s a set.
 */

void pretty_print_set(FILE *fp, set_t s);

END_C_DECL

#endif




/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/





