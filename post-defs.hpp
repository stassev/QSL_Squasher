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


#define  MAX_RES_BITS 20 // Used for generating the Hilbert keys. 
// The largest resolution cube that 
// can be indexed with Hilbert keys has 2^(MAX_RES_BITS*3) elements in 3d.
//**********************************************************************
#define nmax (1024*16) //A dummy set to a value larger than any of NX,NY,NZ.
//**********************************************************************


#define  MAX_RES (1<<MAX_RES_BITS) 	 //It should equal 2^MAX_RES_BITS

#define SIZES()\
    const size_t nx=NX;\
    const size_t ny=NY;\
    const size_t nz=NZ;

#if GEOMETRY==SPHERICAL
	#define solar_radius SOLAR_RADIUS // Don't touch!
	#define solar_radius_sub SOLAR_RADIUS  // Don't touch!
	#define to_radians 1.745329251994329576923690768488612713e-2 // pi/180
#else
	#define solar_radius 1.0 // Don't touch!
    #define solar_radius_sub 0.0 // Don't touch!
    #define to_radians 1.0
#endif

#if QSL_DIM==3
    double LX=(XMAX-XMIN)*to_radians; //Don't touch!
    double LY=(YMAX-YMIN)*to_radians; //Don't touch!
    double ORIGIN_XY[]={XMIN*to_radians,YMIN*to_radians}; //Don't touch!
#else
	#if SLICE_TYPE==SPHERICAL
			double SLICE_NORMAL[] = {0,0,1};   // Don't touch! It's at fixed **radius**
			double SLICE_UP[] = {0,1,0};       // Don't touch!  UP points along latitude
	#endif
#endif

double xmax_file;
double xmin_file;
double ymax_file;
double ymin_file;
double zmax_file;
double zmin_file;

double zmax_fileMid;
double zmin_fileMid;
double ymax_fileMid;
double ymin_fileMid;
double xmax_fileMid;
double xmin_fileMid;



#define INTERP(...) interp_trilinear (__VA_ARGS__)
#define INTERPMID(...) interp_trilinearMid (__VA_ARGS__)
#define INTERP_DIFF(...) interp_trilinear_diff (__VA_ARGS__)



// default integration range for each integration iteration:
#ifndef LOCAL_Q
    #ifndef INTEGRATION_RANGE
        #define INTEGRATION_RANGE 10.0//Mm
    #endif
    #ifndef INTEGRATION_STEPS_PER_CELL
        #define INTEGRATION_STEPS_PER_CELL 5.
    #endif
#endif

bool SAMPLER_INITIALIZED=false;

#if QSL_DIM==2
    double SLICE_RIGHT[] = {0,0,0};//dummy
#endif

#if CALCULATE==QSL
    #define NUM_ODE 10
#else
    #define NUM_ODE 8
#endif



#if QSL_DIM==3
    const size_t init_size=nx_init*ny_init*nz_init;
#endif
#if QSL_DIM==2
    const size_t init_size=nx_init*ny_init;
#endif
