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
 * File: main.cpp
 *       primary command line interface
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Output the title of the software
void titleOutput()
{
		printf("\n");
		printf("###################################################################################\n");
		printf("##  VASP %2.1f Copyright (c) 2014 Brian Y. Chen                                     ##\n", VERS_NUM);
		printf("##                                                                               ##\n");
		printf("##  This program comes with ABSOLUTELY NO WARRANTY. This is free software.       ##\n");
		printf("##  You are welcome to redistribute it under the terms of the GNU General        ##\n");
		printf("##  Public License as published by the Free Software Foundation, either          ##\n");
		printf("##  version 3 of the License, or (at your option) any later version.             ##\n");
		printf("##                                                                               ##\n");
		printf("##  Reference: \"VASP: A Volumetric Analysis of Surface Properties                ##\n");
		printf("##  Yields Insights into Protein-Ligand Binding Specificity.\"                    ##\n");
		printf("##  Brian Y. Chen and Barry Honig, PLOS Computational Biology, 2010.             ##\n");
		printf("###################################################################################\n");
		printf("\n");
}

			 


// This outputs how the code is used.
void defaultOutput()
{
		printf("Example Usage:\n");
		printf("\n");
		printf("PROMPT> ./vasp (-?) \n");
		printf("Returns this basic information. -?, optional, returns verbose information\n");
		printf("\n");
		printf("PROMPT> ./vasp -format\n");
		printf("Generates a description of SURF files used by VASP.\n");
		printf("\n");
		printf("PROMPT> ./vasp -csg <input1.SURF> <input2.SURF> <CSG_OP> <output.SURF> <resolution>\n");
		printf("Perform the CSG operation <CSG_OP> on input1 and input2, at <resolution>\n");
		printf("\n");
		printf("#########################################################################\n");
		printf("\n");
		exit(0);
}

// This shows up when you provide the wrong switch or number of inputs
void errorOutput()
{
		printf("\n");
		printf("===============================================================\n");
		printf("Erroneous Input Detected.  Rerun with no command line arguments\n");
		printf("for usage description and help.\n");
		printf("===============================================================\n");
		printf("\n");
		exit(0);
}

// Detailed usage and explanation
void outputVerboseUsage()
{
		printf("Verbose Usage:\n");
		printf(" The -csg switch:\n");
		printf("     The CSG switch implements three common operations typical of \n");
		printf("     Computational Solid Geometry.  Given two input SURF files, which\n");
		printf("     represent closed surfaces, the -csg switch computes a volumetric\n");
		printf("     intersection (I), union (U), or difference (D) between the\n");
		printf("     inputs, at a given resolution.\n");
		printf("\n");
		printf("     Example> ./vasp -csg 1abc.SURF 1def.SURF D output.SURF .5;\n");
		printf("\n");
		printf(" Typical Usage Recommendations:\n");
		printf("     Recommended resolutions for both switches is .5 Angstroms.  \n");
		printf("     However, finer resolutions can also be used.  Increases in \n");
		printf("     resolution cause increases in memory requirements in an inverse \n");
		printf("     cubic proportion tothe resolution.  Memory usage for the typical\n");
		printf("     protein surface at  .5  Angstroms should peak at several hundred\n");
		printf("     megabytes.  Runtimes at this resolution should be between 15 and\n");
		printf("     60 seconds.\n");
		printf("\n");
		printf(" Reading the text generated in Stdout\n");
		printf("     VASP generates distinctive text output in it's various stages.\n");
		printf("     The nature of the data is described on the left:  either \"STATE\"\n");
		printf("     which indicates what VASP is currently doing, \"INFO\", which\n");
		printf("     provides data about what VASP is currently calculating,  \"ERROR\"\n");
		printf("     which describes what went wrong and in which functions (often\n");
		printf("     followed by exitting the program),  or \"WARNING\" which describes\n");
		printf("     how some kind of contingency analysis was performed because \n");
		printf("     something unusual happened.   Next to the left column, in square\n");
		printf("     brackets, a description informs the user what the STATE, INFO,\n");
		printf("     ERROR, or WARNING is. Finally, on the rest of the line after the\n");
		printf("     brackets, an explanation or the data is provided. Generally, the\n");
		printf("     text generated by VASP doesn't need to be read by users, serving\n");
		printf("     only to inform experts or the curious what the program is doing.\n");
		printf("\n");
		printf(" Referencing VASP\n");
		printf("     VASP was implemented in C++ by Brian Y. Chen, at the Honig Lab,\n");
		printf("     in 2010 and further enhanced by Dr. Chen in his own lab, the \n");
		printf("     Informatics Lab at Lehigh University, through 2014.  If you use\n");
		printf("     VASP to gether results in your research, please cite:\n");
		printf("     \"VASP: A Volumetric Analysis of Surface Properties Yields\n");
		printf("     Insights into Protein-Ligand Binding Specificity\", by Brian Y.\n");
		printf("     Chen and Barry Honig, and also \"VASP-E: Specificity Annotation\n");
		printf("     with a Volumetric Analysis of Electrostatic Isopotentials, by Brian\n");
		printf("     Y. Chen, with which this code release is associated.\n");
		printf("\n");
		printf("     Thank you for doing science.\n");
		printf("#########################################################################\n");
		exit(0);
}

// Detailed explanation of the SURF file format
void outputFileFormat()
{
		printf("SURF FORMAT:\n");
		printf("     -Anywhere in the file, empty lines and line starting with '#' are\n");
		printf("     considered comments, and are ignored.\n");
		printf("     -File has 2 sections: geometry and topology\n");
		printf("     -NO COMMENTS INSIDE EACH SECTION\n");
		printf("\n");
		printf("     The geometry section is always first, and starts with:\n");
		printf("     GEOMETRY: <int>\n");
		printf("     where <int> is the number of points to be specified.\n");
		printf("     then there are <int> lines as follows:\n");
		printf("     <float> <float> <float> <float> <float> <float>\n");
		printf("     which stand for x,y,z, and xnormal, ynormal, znormal.\n");
		printf("\n");
		printf("     The topology section is always second, and starts with:\n");
		printf("     TOPOLOGY: <int>\n");
		printf("     where <int> specifies the number of triangles on the surface\n");
		printf("     Then there are <int> lines as follows:\n");
		printf("     <int> <int> <int>\n");
		printf("     Each int stands for the 3 corners of the triangle\n");
		printf("     and is an index into the array of points provided in the\n");
		printf("     Geometry section.  I.e. the geometry section is indexed\n");
		printf("     starting at zero, and ending at (size-1).\n");
		printf("\n");
		printf("A very simple SURF file example (not including ===):\n");
		printf("=========================================\n");
		printf("####################\n");
		printf("GEOMETRY: 3\n");
		printf("0.0 0.0 0.0 0.0 0.0 1.0\n");
		printf("1.0 0.0 0.0 0.0 0.0 1.0\n");
		printf("0.0 1.0 0.0 0.0 0.0 1.0\n");
		printf("####################\n");
		printf("TOPOLOGY: 1\n");
		printf("0 1 2\n");
		printf("=========================================\n");
		printf("\n");
		printf("#########################################################################\n");
		exit(0);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







//################################################################################################################################
//################################################################################################################################
//################################################################################################################################
//###   #####   #####   #####         ##   #####   ########   #####   ##         ##         ##   ####   ####      ####       #####
//###    ###    ####     ####         ##    ####   ########    ###    ##         ##         ##   ####   ###        ###        ####
//###     #     ###   #   ######   #####     ###   ########     #     ##   ###########   #####   ####   ##    ##    ##   ##    ###
//###           ##   ###   #####   #####      ##   ########           ##       #######   #####          ##   ####   ##   ###   ###
//###   #   #   ##         #####   #####   #   #   ########   #   #   ##       #######   #####          ##   ####   ##   ###   ###
//###   ## ##   ##         #####   #####   ##      ########   ## ##   ##   ###########   #####   ####   ##   ####   ##   ###   ###
//###   #####   ##   ###   #####   #####   ###     ########   #####   ##   ###########   #####   ####   ##    ##    ##   ##    ###
//###   #####   ##   ###   ##         ##   ####    ########   #####   ##         #####   #####   ####   ###        ###        ####
//###   #####   ##   ###   ##         ##   #####   ########   #####   ##         #####   #####   ####   ####      ####       #####
//################################################################################################################################
//################################################################################################################################
//################################################################################################################################
int main(int argc, char* argv[])
{
	
	/////////////////////////////////////////////////
	// Seed the Random Number Generator//////////////
	//
	time_t seconds = time(NULL);
	srandom(seconds);
	
	printf("Random Seed: [%ld]\n", seconds);
	//
	/////////////////////////////////////////////////
	

//	testDebug();


	/////////////////////////////////////////////////
	// Output the title /////////////////////////////
	//
	titleOutput();	
	//
	/////////////////////////////////////////////////

	/////////////////////////////////////////////////
	// default output ///////////////////////////////
	//
	if(argc == 1){ defaultOutput(); }
	//
	/////////////////////////////////////////////////

	/////////////////////////////////////////////////
	// verbose output ///////////////////////////////
	//
	if( argc == 2 && strcmp(argv[1], "-?")==0 ){ outputVerboseUsage(); }
	//
	/////////////////////////////////////////////////

	/////////////////////////////////////////////////
	// file format output ///////////////////////////
	//
	if( argc == 2 && strcmp(argv[1], "-format")==0 ){ outputFileFormat(); }
	//
	/////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// vasp csg //////////////////////////////////////////////////////////////////////////////////////
	// vasp -csg <input1.SURF> <input2.SURF> <CSG_OP> <output.SURF> <resolution>
	//
	if( argc == 7 && strcmp(argv[1], "-csg")==0 ){
		char * input1 = argv[2];
		char * input2 = argv[3];
		char * CSGop = argv[4];
		char * outputSURF = argv[5];
		char * resolution = argv[6];
		
		int op;
		switch(CSGop[0]){
			case 'D':		op = CSG_DIFFERENCE;	break;
			case 'I':		op = CSG_INTERSECTION;	break;
			case 'U':		op = CSG_UNION;		break;
			default :
				printf("ERROR: unrecognized CSG operation: [%s]\n", CSGop);
				exit(1);
		}

		SurfaceObject * surf1 = parseGeometryFile(input1);			// parse the surfaces
		SurfaceObject * surf2 = parseGeometryFile(input2);
		CSG * myCSG = new CSG( surf1, surf2 );						// allocate the CSG object
		SurfaceObject * output = myCSG->runCSG( op, atof(resolution) );	// Run the CSG operation
		generateSURF(output, NULL, outputSURF);						// output the result surface
		
		delete(myCSG);											// clean up
		delete(output);
		delete(surf1);
		delete(surf2);
		exit(1);
	}
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////
	//
	else{ errorOutput(); }
	//
	/////////////////////////////////////////////////
	
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/



