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
 * File: progressBar.cpp
 *       Progress Bar for STDOUT interface
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

//////////////////////////////////////////////////////////////////////
// progressBar.cpp: implementation of Progess Meter
//
//////////////////////////////////////////////////////////////////////

#include "progressBar.h"

void stdoutProgressBar(int curr, int max)
{
        static int progress_count;
        //static int percentage_step = 10;
        static int time_step = 1;
        static char outputString[40];
        static char spinnerString[40];
        static time_t startTime;
        int i = 0;
                
        if(curr==0){ 
                startTime = time(NULL);
                progress_count = 0; 
                spinnerString[0] = '|';
        }
        if(curr == max-1){
                if( curr != 0){ 
                        for(i = 0; i<((int)strlen(outputString)); i++){ printf("\b");   }
                }                       
                sprintf(outputString, "Progress: 100.0%%\n");
                printf("%s", outputString);
                fflush(stdout);
        }

        ///progress time based
        time_t currentTime = time(NULL);
        if( currentTime-startTime >= progress_count ){
                if( curr != 0){ 
                        for(i = 0; i<((int)strlen(outputString)); i++){ printf("\b");   }
                }
                switch(spinnerString[0]){
                        case '|': spinnerString[0] = '/'; break;
                        case '/': spinnerString[0] = '-'; break;
                        case '-': spinnerString[0] = '\\'; break;
                        case '\\': spinnerString[0] = '|'; break;
                }
                sprintf(outputString, "Progress: %5.1f%% %s   ", ((double)(100*curr))/((double)max), spinnerString );
                printf("%s", outputString);
                fflush(stdout);
                progress_count += time_step;
                startTime = currentTime;
        }

        ///progress increment based
        /*
        if( ((double)(100*curr))/((double)max) >= progress_count ){
                if( curr != 0){ 
                        for(i = 0; i<strlen(outputString); i++){ printf("\b");   }
                }                       
                sprintf(outputString, "Progress: %5.1lf%%", ((double)(100*curr))/((double)max) );
                printf("%s", outputString);
                fflush(stdout);
                progress_count += percentage_step;
        }
        */

        
}





/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/





