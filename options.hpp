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

#include "pre-defs.hpp"

//#define OpenCL_DEVICE_TYPE CL_DEVICE_TYPE_GPU
#define OpenCL_DEVICE_TYPE CL_DEVICE_TYPE_CPU


#define SOLAR_RADIUS 696.0

//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************


 
std::string in_dir  ="./cartesian_demo/";
std::string in_filename="";
 
#define NX 128
#define NY 128
#define NZ 128
 
//// Your coordinates are (pick one): 
#define GEOMETRY CARTESIAN
//#define GEOMETRY SPHERICAL
//**********************************


#if GEOMETRY==CARTESIAN
	#define SLICE_TYPE CARTESIAN // Don't change!
#else
	//#define SLICE_TYPE CARTESIAN
	#define SLICE_TYPE SPHERICAL
#endif


//**********************************************************************

#define NGPU 0

#define INTEGRATION_SCHEME EULER 
//#define INTEGRATION_SCHEME ADAPTIVE
#if INTEGRATION_SCHEME==ADAPTIVE
	#define eps_rel 1.0e-2
	#define eps_abs 1.0e-2
#endif 
//**********************************************************************

#define MAX_HALF_LENGTH_FIELD_LINE (SOLAR_RADIUS*10.) //Mm 

// #define LOCAL_Q
// #define INTEGRATION_RANGE 10.
#define INTEGRATION_STEPS_PER_CELL 5.

#ifndef LOCAL_Q
	#define MARK_OPEN_FIELD_LINES 
#endif

//#define CALCULATE FIELD_LINE_LENGTH
#define CALCULATE QSL

//**********************************************************************
//**********************************************************************
//**********************************************************************
// Do you want a 2d slice (QSL_DIM=2) of Q values, or do you prefer to 
// calculate the Q values in 3d (QSL_DIM=3)?
//#define QSL_DIM 3
#define QSL_DIM 3
//**********************************************************************
// Initial grid size for Q before refining along Hilbert curve.
#if QSL_DIM==2 
	const size_t nx_init = 1024;
	const size_t ny_init = 1024;
#endif
#if QSL_DIM==3
	const size_t nx_init = 64;
	const size_t ny_init = 64;
	const size_t nz_init = 64;
#endif

#if QSL_DIM==2
	#define ZMIN  0.001    //Mm
    #if SLICE_TYPE==CARTESIAN
		double SLICE_NORMAL[] = {0,1,0};   
		double SLICE_UP[] = {0,0,1}; 
		double SLICE_CENTER[] = {0,0,7.15}; // Mm
		double SLICE_LX = 18.; //Mm
		double SLICE_LY = 14.; //Mm
	#endif
    #if SLICE_TYPE==SPHERICAL // SPHERICAL slice at fixed r
		double SLICE_CENTER[] = {0,0,10.}; // {lon in deg, lat in deg, r in Mm}
		double SLICE_LX = 18.;  //range in longitude in deg
		double SLICE_LY = 28.;  //range in latitude in deg
	#endif
#endif
#if QSL_DIM==3
    #if GEOMETRY==CARTESIAN
		#define XMIN -9.9      //Mm
        #define XMAX  9.9      //Mm
        #define YMIN -14.9     //Mm
        #define YMAX  14.9     //Mm
        #define ZMIN  0.00001  //Mm
        #define ZMAX  14.3     //Mm
        #define z_sampler(z) \
			z=z*(ZMAX-ZMIN)+ZMIN;
    #else
		#define XMIN -8.       //degrees
		#define XMAX  8.       //degrees
		#define YMIN -12.      //degrees
		#define YMAX  12.      //degrees
		#define ZMIN  0.00001  //Mm
		#define ZMAX  143.     //Mm
		#define z_sampler(z) \
			z=(exp(pow((70.+z*(300.))/370.,3))-exp(pow(70./370.,3)))/(exp(1.)-exp(pow(70./370.,3)))*(ZMAX-ZMIN)+ZMIN;
	#endif
#endif



// Pick one(!) interpolation algorithm used for interpolating 
// the B-field values:
#define INTERPOLATION_TYPE TRILINEAR
//#define INTERPOLATION_TYPE TRIQUADRATIC
//#define INTERPOLATION_TYPE TRICUBIC
//**********************************************************************

//If jump is in the log10(length), then a good guess here is 0.05
#define LENGTH_JUMP_REFINEMENT_THRESHOLD 1.0 //Mm

#if INTEGRATION_SCHEME==EULER
	const size_t CHUNKSIZE = pow(2,19);
#else
    const size_t CHUNKSIZE = pow(2,16);
#endif
//**********************************************************************
#define  MAX_RES_BITS 20 // Used for generating the Hilbert keys. 
// The largest resolution cube that 
// can be indexed with Hilbert keys has 2^(MAX_RES_BITS*3) elements in 3d.
//**********************************************************************
#define nmax (1024*16) //A dummy set to a value larger than any of NX,NY,NZ.
//**********************************************************************

#define DISPLACEMENT_WEIGHT 10. // This factor gives more weight to derivatives when calculating errors for adaptive stepper



#include "post-defs.hpp"
