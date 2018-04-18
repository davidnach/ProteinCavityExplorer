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
 * File: flunclib.h
 *       implementation of minor functions
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

#include "StdAfx.h"



/////////////////////////////
///We sort objs, and re-arrange refs in the same way we re-arranged objs.
///ascending
void smallSort(int * objs, int * refs, int size){
	int min = INT_MAX;
	int tempRef = 0;
	int tempVal;
	int i = 0;
	int j = 0;

	for(i = 0; i<size; i++){
		//Reset the min val, to keep finding the min
		min = INT_MAX;
		//find the smallest value between the val after i to the end
		for(j = i; j<size; j++){
			if(objs[j] <= min){ 
				tempRef = j; 
				min = objs[j];
			}
		}
		//reArrange the objs array;
		tempVal = objs[i];
		objs[i] = objs[tempRef];
		objs[tempRef] = tempVal;
		//reArrange the refs array;
		tempVal = refs[i];
		refs[i] = refs[tempRef];
		refs[tempRef] = tempVal;
	}
}


/////////////////////////////
///We sort objs, and re-arrange refs in the same way we re-arranged objs.
///ascending
void smallSort(double * objs, int * refs, int size){
	double min = HUGE_VAL;
	int tempRef = 0;
	double tempVal1;
	int tempVal2;
	int i = 0;
	int j = 0;

	for(i = 0; i<size; i++){
		//Reset the min val, to keep finding the min
		min = HUGE_VAL;
		//find the smallest value between the val after i to the end
		for(j = i; j<size; j++){
			if(objs[j] <= min){ 
				tempRef = j; 
				min = objs[j];
			}
		}
		//reArrange the objs array;
		tempVal1 = objs[i];
		objs[i] = objs[tempRef];
		objs[tempRef] = tempVal1;
		//reArrange the refs array;
		tempVal2 = refs[i];
		refs[i] = refs[tempRef];
		refs[tempRef] = tempVal2;
	}
}


/////////////////////////////
///We sort objs, and re-arrange refs in the same way we re-arranged objs.
void smallSortDescending(int * objs, int * refs, int size){
	int max = 0;
	int tempRef = 0;
	int tempVal;
	int i = 0;
	int j = 0;

	for(i = 0; i<size; i++){
		//Reset the min val, to keep finding the min
		max = 0;
		//find the smallest value between the val after i to the end
		for(j = i; j<size; j++){
			if(objs[j] >= max){ 
				tempRef = j; 
				max = objs[j];
			}
		}
		//reArrange the objs array;
		tempVal = objs[i];
		objs[i] = objs[tempRef];
		objs[tempRef] = tempVal;
		//reArrange the refs array;
		tempVal = refs[i];
		refs[i] = refs[tempRef];
		refs[tempRef] = tempVal;
	}
}


///ascending
void mergeSort(int * objs, int * refs, int size)
{
	int i = 0;
	double * sorter = new double[size];
	for(i = 0; i<size; i++){
		sorter[i] = (double) objs[i];
	}

	mergeSort(sorter, refs, size);
	for(i = 0; i<size; i++){
		objs[i] = (int) sorter[i];
	}

	delete[](sorter);
}

///ascending
void mergeSort(double * objs, int * refs, int size)
{
	///if its big, split.
	if(size > 10){
		int size1 = ((int) ( ((double)size)/2.0) );
		int size2 = size - size1;

		double * objs1 = new double[size1];
		int * refs1 = new int[size1];
		double * objs2 = new double[size2];
		int * refs2 = new int[size2];

		int i = 0;

		////split into two problems.
		for(i = 0; i<size1; i++){
			objs1[i] = objs[i];
			refs1[i] = refs[i];
		}
		for(i = 0; i<size2; i++){
			objs2[i] = objs[size1+i];
			refs2[i] = refs[size1+i];
		}

		///sort them
		mergeSort(objs1, refs1, size1);
		mergeSort(objs2, refs2, size2);

		///merge the two solutions
		int counter1 = 0;
		int counter2 = 0;
		int caseFlag;
		for(i = 0; i<size; i++){
			caseFlag = 0;
			if(counter1 < size1 && counter2 < size2){
				caseFlag = 0;
			}else if(counter1 < size1 && counter2 >= size2){
				caseFlag = 1;
			}else if(counter1 >= size1 && counter2 < size2){
				caseFlag = 2;
			}else {
				caseFlag = 3;	///error
			}
			switch(caseFlag){
				case 0 : 
					if( (objs1[counter1] <= objs2[counter2]) ){
						objs[i] = objs1[counter1];
						refs[i] = refs1[counter1];
						counter1++;
					}else{
						objs[i] = objs2[counter2];
						refs[i] = refs2[counter2];
						counter2++;
					}break;
				case 1 :
					objs[i] = objs1[counter1];
					refs[i] = refs1[counter1];
					counter1++;
					break;
				case 2 :
					objs[i] = objs2[counter2];
					refs[i] = refs2[counter2];
					counter2++;
					break;
				case 3 :
					printf("ERROR: mergeSort is buggy\n");
					exit(0);
			}
		}

		delete[](objs1);
		delete[](refs1);
		delete[](objs2);
		delete[](refs2);

	}
	else{////otherwise insertion sort it.
		smallSort(objs, refs, size);		
	}




}








