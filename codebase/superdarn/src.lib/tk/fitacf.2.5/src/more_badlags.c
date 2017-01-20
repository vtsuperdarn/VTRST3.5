/* more_badlags.c
   ==============
   Author: R.J.Barnes & K.Baker & P.Ponomarenko
*/

/*
 LICENSE AND DISCLAIMER

 Copyright (c) 2012 The Johns Hopkins University/Applied Physics Laboratory

 This file is part of the Radar Software Toolkit (RST).

 RST is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 RST is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with RST.  If not, see <http://www.gnu.org/licenses/>.



*/



#include <math.h>
#include <stdio.h>

double more_badlags(double *w,int *badlag,
                   double noise_lev,int mplgs,int nave) {

    double fluct0, fluct, fluct_old;
    short int badflag_1, badflag_2, k_old, k;
    short int sum_np;

    badflag_1 = 0;
    badflag_2 = 0;

    fprintf(stderr,"In more_badlags.c\n");
    fprintf(stderr,"sum_np starts as:    %d\n", sum_np);
    fprintf(stderr,"mplgs passed in as:  %d\n", mplgs);

    fluct0 =  w[0]/sqrt(2.0*(double) nave);
    fluct =  w[0] + 2.0*noise_lev+fluct0;
    fluct_old = fluct;
    sum_np = 0;
    k_old = 0;

    for (k=0; k<mplgs; k++) {
/*        fprintf(stderr, "top badlag[k]:  %d\n",badlag[k]); */
        if (badlag[k]){
            fprintf(stderr, "Found badlag, continuing on now.  k: %d\n", k);
            continue;
        }
        fprintf(stderr,"More badlags loop step 1\n");
        if (badflag_2) badlag[k]=7;
        else if (w[k] <= w[0]/sqrt((double) nave)) {  /* if (w[k] <= 1.0) { */
            fprintf(stderr,"More badlags loop step 2\n");
   	        badlag[k] = 3;
   	        badflag_2 = badflag_1;
   	        badflag_1 = 1;
        } else {
            fprintf(stderr, "more badlags loop step 3\n");
            badflag_1 = 0;
   	        if (w[k] > fluct) {
	            badlag[k] = 5;
	            if (k < (mplgs - 1)) {
		            if ((w[k] < fluct_old) && (w[k+1] > fluct) && (w[k+1] < w[k])) {
		                badlag[k_old] = 9;
		                --sum_np;
		                badlag[k] = 0;
		            }
                }
	        }
            fluct_old = fluct;
            fluct = 2.0*noise_lev + w[k] + fluct0;
        }
        fprintf("badlag[k]:    %d\n", badlag[k]);
        if (!badlag[k]) {
            ++sum_np;
	        k_old = k;
        }
    }
    fprintf(stderr,"sum_np in more badlags: %d\n",sum_np);
    return (double) sum_np;
}
