/****************************************************************************
 * SurfProcessing
 *   Copyright (c) 2014 Brian Y. Chen
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
 * This file incorporates work covered by the following copyright
 * and permission notices:
 *
 * Copyright (c) 1998-2014 Geometric Tools, LLC
 *
 * Permission to use, copy, modify, and/or distribute this software
 * for any purpose with or without fee is hereby granted, provided
 * that the above copyright notice and this permission notice appear
 * in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
 * AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFTS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * File: triIntersect.cpp
 *       interface for geometric tests for triangles
 *
 * Written by 
 *       David Eberly
 * Modified by 
 *       Brian Y. Chen
 *
 ***************************************************************************/

#ifndef POINT_TRI_DIST_H
#define POINT_TRI_DIST_H

#define VEC( v, p1, p2) \
	v[0] = p2[0]-p1[0];	\
	v[1] = p2[1]-p1[1];	\
	v[2] = p2[2]-p1[2];	\



#include "StdAfx.h"

///this returns the actual distance
//double distanceToSegment(double * A, double * B, double * P);
///this returns the square of the distance
double distanceToSegmentFast(double * A, double * B, double * P);












#endif



