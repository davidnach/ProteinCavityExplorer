==============================================================================
This document describes data and source code used as input for the following 
publication:

VASP-E: Specificity Annotation with a Volumetric Analysis of
Electrostatic Isopotentials

by Brian Y. Chen, 
Dept. Computer Science and Engineering, Lehigh University

Please direct any questions to <chen@cse.lehigh.edu>
==============================================================================

License:

The source code provided here is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

The source code provided is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
for more details.

A copy of the GNU General Public License is available in each software 
directory described below.  If you have not received this license, see 
<http://www.gnu.org/licenses/>.  


Compilation:

The software below are provided in two folders: ./surfProcessing/ and ./vasp/
To compile the code, change your working directory to one of these folders 
and type "make".  The executable binary will be generated in the "debug" folder.
Object files are generated in the "object" folder.  Source code is in the
"source" folder.



File Format:

In addition to PDB file formats, the software provided also operates on closed
surfaces of triangles described by .SURF files.  The .SURF file format is very
simple, and a concise description of the format is provided in here:

SURFformat.txt


Capabilities of vasp:

This program, an acronym for "Volumetric Analysis of Surface Properties", is
designed to perform Constructive Solid Geometry operations based on .SURF
files.

./vasp -csg <input1.SURF> <input2.SURF> <CSG_OP> <output.SURF> <resolution>

It accepts two input surfaces in the .SURF format, substituting for
<input1.SURF> and <input2.SURF>, and performs one of three CSG operations: 
Intersection ("I"), Union ("U"), or Difference ("D").  The output is
generated in a .SURF file named <output.SURF>, even if the output has no
surface at all (for example, the intersection of two objects that do not
overlap).  The result of the final surface is estimated at <resolution>,
a value for which we recommend .5 angstroms for performance reasons on
modern computers with standard configurations.

Smaller values of <resolution> will produce finer, more accurate surfaces
at the expense of longer computation times and larger memory usage.

The executable ./vasp may be used to operate on any .SURF file, whether it
is was generated as a molecular surface or as a electrostatic isopotential.

***************************************************************************
If you generate your own .SURF files it is crucial to note that the surface
must be mathematically closed, otherwise the interior and exterior of the
surface are undefined.  Even the smallest openings can create errors.  If 
a nonclosed surface is provided as input, VASP may crash ungracefully, 
reporting the failure of certain geometric tests.
***************************************************************************


Unit Tests:

Five surfaces are provided for unit testing.  These surfaces are:

test1.SURF        # a sphere of radius 3, centered at the origin
test2.SURF        # a sphere of radius 3, centered at 1,0,0
difference.SURF   # the CSG difference of test1.SURF and test2.SURF
union.SURF        # the CSG union of test1.SURF and test2.SURF
intersection.SURF # the CSG intersection of test1.SURF and test2.SURF

You can visualize these surfaces, and any .SURF file, using SURFview, 
which is part of this release.

test1 and test2 are spheres that are offset slightly; when viewed together
they look like a typical Venn diagram.  The difference, union, and
intersection surfaces should reflect the expected appearance, except that
the difference may have a jagged, noisy edge, especially if computed at .5 
angstroms resolution.  This effect occurs because of the very thin surface
generated at the edge of the hollowed out sphere.  For molecular surfaces,
and binding cavities, which are considerably larger, these noisy effects
have a trivial effect on accuracy.  Alternatively, you can change the 
resolution parameter to a smaller threshold, such as .05 to create a smooth
edge for these small surfaces.  Such fine resolutions are not recommend
for larger surfaces, such as molecular surfaces and cavities, because they
will require too much memory.

If VASP reports a segfault immediately upon running, you may be requesting
too fine a resolution for your surface.  Use a coarser surface.






