/*
	This file is part of QSL Squasher. 
	Copyright (C) 2014-2019  Svetlin Tassev
							 Harvard-Smithsonian Center for Astrophysics
							 Braintree High School
	
     QSL Squasher is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
*/

VEX_FUNCTION(cl_ushort4, find_index, (double, x)(double, y)(double,z)(cl_double4*, ff)(cl_ushort4*, ff_lookup),
    /*
     * Finds the index in the B arrays where the point {x,y,z} is located.
     */
    
    
    DEFS();
    
    double xmin=ff[0].s0;
	double xmax=ff[nx-1].s0;
	double ymin=ff[0].s1;
	double ymax=ff[ny-1].s1;
	double zmin=ff[0].s2;
	double zmax=ff[nz-1].s2;
	
	// Check for out-of-range values
    ushort a=1;
    
	#ifdef PERIODIC_XY
		x=fmod(x-ff[0].s0+ff[nx-1].s0-ff[0].s0,ff[nx-1].s0-ff[0].s0)+ff[0].s0;
		y=fmod(y-ff[0].s1+ff[ny-1].s1-ff[0].s1,ff[ny-1].s1-ff[0].s1)+ff[0].s1;
	#else
		#ifndef GLOBAL_MODEL
		if (x<=xmin) {
			x=xmin+1.e-4;
				a=0;
		}
		if (x>=xmax)  {
			x=xmax-1.e-4;
			a=0;
		}
		#else
			x=fmod((x+M_PI*2.0000001-xmin),(2.*M_PI))+xmin;
		#endif
		
		#ifdef GLOBAL_MODEL
			if (fabs(y)>1.5706) 
				a=0;
		#endif
		
		if (y<=ymin) {
			y=ymin+1.e-4; 
			#ifndef GLOBAL_MODEL
				a=0;
			#endif
		}
		if (y>=ymax) {
			y=ymax-1.e-4;
			#ifndef GLOBAL_MODEL
				a=0;
			#endif
		}
	#endif


    
    if (z<=z_minimum) {
        z=z_minimum+1.e-4;
        a=0;
    }

    if (z>=zmax)  {
         z=zmax-1.e-4;
         a=0;
    }
    
    ushort i;
    ushort j;
    ushort k;
    
    //Guess the index by using lookup tables
    //i =(x-ff[0].s0)/(ff[1].s0-ff[0].s0)+0.5; 
    //j =(y-ff[0].s1)/(ff[1].s1-ff[0].s1)+0.5;
    //k =(z-ff[0].s2)/(ff[1].s2-ff[0].s2)+0.5;
    
    ushort ix=((x-xmin)/(xmax-xmin))*((double)nmax)+0.5;
    ushort iy=((y-ymin)/(ymax-ymin))*((double)nmax)+0.5;
    ushort iz=((z-zmin)/(zmax-zmin))*((double)nmax)+0.5;
    
    if (ix>nmax) ix=nmax-1;
    if (iy>nmax) iy=nmax-1;
    if (iz>nmax) iz=nmax-1;
    
    i=ff_lookup[ix].s0;
    j=ff_lookup[iy].s1;
    k=ff_lookup[iz].s2;
    
    //Don't want to overshoot:
    if (i>nx-2) i=nx-2;
    if (j>ny-2) j=ny-2;
    if (k>nz-2) k=nz-2;
   
    //Starting with the guess above, find the exact index location
    #ifdef GLOBAL_MODEL
		while ((i<=nx-1) && (ff[i].s0<x))i++;
    #else
		while ((i<nx-1) && (ff[i].s0<x))i++;
    #endif
    while ((j<ny-1) && (ff[j].s1<y))j++;
    while ((k<nz-1) && (ff[k].s2<z))k++;
    
    while ((i>=0) && (ff[i].s0>x))i--;
    while ((j>=0) && (ff[j].s1>y))j--;
    while ((k>=0) && (ff[k].s2>z))k--;
    
    
    
    #ifdef GLOBAL_MODEL
		if (x>=ff[nx-1].s0) i=nx-1;
    #endif
    
    ushort4 ind={i,j,k,a};
    return ind;
);



#if CALCULATE!=QSL
VEX_FUNCTION(cl_ushort4, find_indexMid, (double, x)(double, y)(double,z)(cl_double4*, ff)(cl_ushort4*, ff_lookup),
    /*
     * Finds the index in the B arrays where the point {x,y,z} is located.
     */
    
    
    DEFS();
    
    double xmin=ff[0].s0;
	double xmax=ff[nx-2].s0;
	double ymin=ff[0].s1;
	double ymax=ff[ny-2].s1;
	double zmin=ff[0].s2;
	double zmax=ff[nz-2].s2;
	
	// Check for out-of-range values
    ushort a=1;
    
	#ifdef PERIODIC_XY
		x=fmod(x-ff[0].s0+ff[nx-2].s0-ff[0].s0,ff[nx-2].s0-ff[0].s0)+ff[0].s0;
		y=fmod(y-ff[0].s1+ff[ny-2].s1-ff[0].s1,ff[ny-2].s1-ff[0].s1)+ff[0].s1;
	#else
		#ifndef GLOBAL_MODEL
		if (x<=xmin) {
			x=xmin+1.e-4;
				a=0;
		}
		if (x>=xmax)  {
			x=xmax-1.e-4;
			a=0;
		}
		#else
			x=fmod((x+M_PI*2.0000001-xmin),(2.*M_PI))+xmin;
		#endif
		
		#ifdef GLOBAL_MODEL
			if (fabs(y)>1.5706) 
				a=0;
		#endif
		
		if (y<=ymin) {
			y=ymin+1.e-4; 
			#ifndef GLOBAL_MODEL
				a=0;
			#endif
		}
		if (y>=ymax) {
			y=ymax-1.e-4;
			#ifndef GLOBAL_MODEL
				a=0;
			#endif
		}
	#endif


    
    if (z<=z_minimum) {
        z=z_minimum+1.e-4;
        a=0;
    }
    if (z<=zmin) {
        z=zmin+1.e-4;
        a=0;
    }
    if (z>=zmax)  {
         z=zmax-1.e-4;
         a=0;
    }
    
    ushort i;
    ushort j;
    ushort k;
    
    //Guess the index by using lookup tables
    //i =(x-ff[0].s0)/(ff[1].s0-ff[0].s0)+0.5; 
    //j =(y-ff[0].s1)/(ff[1].s1-ff[0].s1)+0.5;
    //k =(z-ff[0].s2)/(ff[1].s2-ff[0].s2)+0.5;
    
    ushort ix=((x-xmin)/(xmax-xmin))*((double)nmax)+0.5;
    ushort iy=((y-ymin)/(ymax-ymin))*((double)nmax)+0.5;
    ushort iz=((z-zmin)/(zmax-zmin))*((double)nmax)+0.5;
    
    if (ix>nmax) ix=nmax-1;
    if (iy>nmax) iy=nmax-1;
    if (iz>nmax) iz=nmax-1;
    
    i=ff_lookup[ix].s0;
    j=ff_lookup[iy].s1;
    k=ff_lookup[iz].s2;
    
    //Don't want to overshoot:
    if (i>nx-3) i=nx-3;
    if (j>ny-3) j=ny-3;
    if (k>nz-3) k=nz-3;
   
    //Starting with the guess above, find the exact index location
    #ifdef GLOBAL_MODEL
		while ((i<=nx-2) && (ff[i].s0<x))i++;
    #else
		while ((i<nx-2) && (ff[i].s0<x))i++;
    #endif
    while ((j<ny-2) && (ff[j].s1<y))j++;
    while ((k<nz-2) && (ff[k].s2<z))k++;
    
    while ((i>=0) && (ff[i].s0>x))i--;
    while ((j>=0) && (ff[j].s1>y))j--;
    while ((k>=0) && (ff[k].s2>z))k--;
    
    #ifdef GLOBAL_MODEL
		if (x>=ff[nx-2].s0) i=nx-2;
    #endif
    
    ushort4 ind={i,j,k,a};
    return ind;
);
#endif
