##
## VASP: Volumetric Analysis of Surface Properties
## Copyright (c) 2014 Brian Y. Chen
## All rights reserved.
##
## This file is part of VASP
## 
## VASP is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## VASP is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with VASP.  If not, see <http://www.gnu.org/licenses/>.  
##
## File: delphi_customCRG.pl
##       A frontend for DelPhi that accepts custom charge values
##
## Written by 
##       Brian Y. Chen <chen@lehigh.edu>
##
## WWW URL: http://cse.lehigh.edu/~chen/
## Email: chen@lehigh.edu
## Documentation can be found here: 
## http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
##
##

#!/usr/bin/perl -w
use strict;

#DelPhi parameters
my $perfil = 25;
my $gsize = 145;
my $bndcon_init = 4;
my $salt = 0.145;
my $prbrad = 1.4;
my $linit = 1000;
my $nonit = 1000;
#Number of focussing steps:
my $num_steps = 3;

#DelPhi executable
my $DELPHI = '/infolab/exchange/software/delphi_new/executable/delphi';

#atom partial charges
my $CRG = $ARGV[1];

#atom radii
my $SIZ = '/infolab/exchange/software/delphi_new/parameters/amber.siz';

my $phi_prev;

for my $step ( 0 .. $num_steps ) {

    print "Commencing focusing step $step of $num_steps ... ";

    my $phi = 'map_perfil' . $perfil . '.phi';
    my $out = 'out_perfil' . $perfil . '.log';

    open ( DELPHI, ">param" )
        or die;

    print DELPHI "gsize=$gsize\n";
    print DELPHI "perfil=$perfil\n";
    print DELPHI "bndcon=$bndcon_init\n" if $step == 0;
    print DELPHI "bndcon=3\n" if $step > 0;
    print DELPHI "salt=$salt\n";
    print DELPHI "prbrad=$prbrad\n";
    print DELPHI "linit=$linit\n" if ( $linit );
    print DELPHI "nonit=$nonit\n" if ( $nonit );
    print DELPHI "in(pdb,file=$ARGV[0])\n";
    print DELPHI "in(crg,file=$CRG)\n";
    print DELPHI "in(siz,file=$SIZ)\n";
    print DELPHI "in(phi,file=$phi_prev)\n" if ( $step > 0 );
    print DELPHI "out(phi,file=" . $phi . ")\n" if ( $step < 3 );
    print DELPHI "energy(s,c,g)\n";
    print DELPHI "out(modpdb,unit=19)\n";
    print DELPHI "out(phi,file=out.grd,format=BIOSYM)\n" if ( $step == 3 );
    print DELPHI "out(modpdb, file=file.modpdb)";  ##newshit

    close DELPHI;

    `$DELPHI param > $out`;

    `cp $out final.log`
        if ( $step == 3 );


    $phi_prev = $phi;

    $perfil += 20;

    print "done.\n";

}


