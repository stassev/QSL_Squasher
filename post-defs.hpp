/*
	This file is part of QSL Squasher. 
	Copyright (C) 2014, 2015, 2016  Svetlin Tassev
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


#if INTERPOLATION_TYPE==TRILINEAR
        #define INTERP(...) interp_trilinear (__VA_ARGS__)
        #define INTERP_DIFF(...) interp_trilinear_diff (__VA_ARGS__)
#endif
#if INTERPOLATION_TYPE==TRIQUADRATIC
        #define INTERP(...) interp_triquadratic (__VA_ARGS__)
        #define INTERP_DIFF(...) interp_triquadratic_diff (__VA_ARGS__)
#endif
#if INTERPOLATION_TYPE==TRICUBIC
        #define INTERP(...)  interp_tricubic (__VA_ARGS__)
        #define INTERP_DIFF(...) interp_tricubic_diff (__VA_ARGS__)
#endif


// default integration range for each integration iteration:
#ifndef LOCAL_Q
    #ifndef INTEGRATION_RANGE
        #define INTEGRATION_RANGE 10.0//Mm
    #endif
    #ifndef INTEGRATION_STEPS_PER_CELL
        #if INTEGRATION_SCHEME==EULER
			#define INTEGRATION_STEPS_PER_CELL 30.
		#else
			#define INTEGRATION_STEPS_PER_CELL 1.
		#endif
    #endif
#endif

bool SAMPLER_INITIALIZED=false;

#if QSL_DIM==2
    double SLICE_RIGHT[] = {0,0,0};//dummy
#endif

#if CALCULATE==QSL
    #define NUM_ODE 10
#else
    #define NUM_ODE 4
#endif

