##***************************************************************************
## SurfProcessing
## Copyright (c) 2014 Brian Y. Chen
## All rights reserved.
##
## This file is part of SurfProcessing
## 
## SurfProcessing is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## SurfProcessing is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SurfProcessing.  If not, see <http://www.gnu.org/licenses/>.  
##
## File: Makefile
##       compilation script
##
## Written by 
##       Brian Y. Chen <chen@lehigh.edu>
##
## WWW URL: http://cse.lehigh.edu/~chen/
## Email: chen@lehigh.edu
## Documentation can be found here: 
## http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
##
##**************************************************************************/



###################################################################################
## ### ##   ## ### #    #    #   # ####    ########################################
##  #  # ### # ## ## #### ##### ## #### ###########################################
## # # # ### # # ### #### ##### ## #### ###########################################
## ### #     #   ###   ##   ### ## ####   #########################################
## ### # ### # ## ## #### ##### ## #### ###########################################
## ### # ### # ### # #### ##### ## #### ###########################################
## ### # ### # ### #    # ####   #    #    ########################################
###################################################################################



###################################################################################
###   # ####  ###   ##  ###########################################################
### ### ### ## # #### #############################################################
###  ## ###    # #  ## ############################################################
### ### ### ## # ## ### ###########################################################
### ###   # ## ##  ##  ############################################################
###################################################################################

###################
#compiler options #
###################
CC = gcc
CXX = g++
CPPFLAGS = -g -pedantic -pedantic-errors -Wall -Wno-deprecated -Wparentheses 
CFLAGS = -g 

###################
#directory settings
###################
SOURCEDIR = source/
OBJECTDIR = object/
DEBUGDIR  = debug/
EXENAME   = surfProcessing
EXEPATH   = $(DEBUGDIR)$(EXENAME)
LIBS	  = 
DEFINES   = 
INCLUDES  = -I./libs/trollbaseHeaders/ \
#\-I./libs/qhullHeaders/



###################################################################################
###   #   # ###   ##  #############################################################
### #### ## ### ### ###############################################################
###  ### ## ###  ### ##############################################################
### #### ## ### ##### #############################################################
### ###   #   #   #  ##############################################################
###################################################################################

C_OBJS = \
$(OBJECTDIR)set.o \
$(OBJECTDIR)prime.o \

BATCH_OBJS = \
$(OBJECTDIR)SurfaceObject.o \
$(OBJECTDIR)mathlib.o \
$(OBJECTDIR)triIntersect.o \
$(OBJECTDIR)primativesProcessing.o \
$(OBJECTDIR)CubeTable.o \
$(OBJECTDIR)LatticeHash.o \
$(OBJECTDIR)SurfaceReconstruction.o \
$(OBJECTDIR)LatticeObj.o \
$(OBJECTDIR)TriangleLatticeHash.o \
$(OBJECTDIR)vertTriSet.o \
$(OBJECTDIR)SurfaceOutput.o \
$(OBJECTDIR)funclib.o \
$(OBJECTDIR)pointTriDist.o \
$(OBJECTDIR)GeometryParser.o \
$(OBJECTDIR)main.o \
$(OBJECTDIR)phimapGrid.o \



###################################################################################
####  ##  ## ### #   ##   # ####  ##   #   ##  ## ## ##############################
### ### ## #  #  # ## ## ## ### ## ## ### ## ## #  # ##############################
### ### ## # # # #   ### ## ###    ## ### ## ## # #  ##############################
### ### ## # ### # ##### ## ### ## ## ### ## ## # ## ##############################
####  ##  ## ### # ####   #   # ## ## ##   ##  ## ## ##############################
###################################################################################

all: $(BATCH_OBJS) $(C_OBJS)
	$(CXX) $(CPPFLAGS) -o $(EXEPATH) $(BATCH_OBJS) $(C_OBJS) $(LIBS) $(INCLUDES)



###################################################################################
# Dependencies                                                                    #
###################################################################################

$(BATCH_OBJS): $(OBJECTDIR)%.o: $(SOURCEDIR)%.cpp
	$(CXX) -c $(CPPFLAGS) $< -o $@ $(INCLUDES)

$(C_OBJS): $(OBJECTDIR)%.o: $(SOURCEDIR)%.c
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDES)



###################################################################################
####  # ###   ##  ## ## ###########################################################
### ### ### ### ## #  # ###########################################################
### ### ###  ##    # #  ###########################################################
### ### ### ### ## # ## ###########################################################
####  #   #   # ## # ## ###########################################################
###################################################################################

clean:
	rm -rf $(BATCH_OBJS) $(C_OBJS) $(EXEPATH) object/*




