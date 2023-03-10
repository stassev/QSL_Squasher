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
VEX_FUNCTION(cl_double4, interp_trilinear, (double, x)(double, y)(double,z)(cl_double4*, B)(cl_double4*, ff)(cl_ushort4, ind),
    DEFS();
    double a=1.0;
    if (ind.s3 == 0) a=1e-9;
    
    #ifdef PERIODIC_XY
		x=fmod(x-ff[0].s0+ff[nx-1].s0-ff[0].s0,ff[nx-1].s0-ff[0].s0)+ff[0].s0;
		y=fmod(y-ff[0].s1+ff[ny-1].s1-ff[0].s1,ff[ny-1].s1-ff[0].s1)+ff[0].s1;
	#else 
		#ifndef GLOBAL_MODEL
			if (x<=ff[0].s0) 
				x=ff[0].s0+1.e-4;
			
			if (x>=ff[nx-1].s0)  
				x=ff[nx-1].s0-1.e-4; 
		#else
			x=fmod((x+M_PI*2.0000001-ff[0].s0),(2.*M_PI))+ff[0].s0;
		#endif
		if (y<=ff[0].s1) 
			y=ff[0].s1+1.e-4; 
		if (y>=ff[ny-1].s1) 
			y=ff[ny-1].s1-1.e-4; 
	#endif
	
    if (z<=z_minimum) 
        z=z_minimum+1.e-4; 
    
    if (z>=ff[nz-1].s2)  
        z=ff[nz-1].s2-1.e-4; 
    
    
    size_t i = ind.s0;
    size_t j = ind.s1;
    size_t k = ind.s2;
    
    double fx=ff[i].s0;
    double fy=ff[j].s1;
    double fz=ff[k].s2;
    
    double xP;
    size_t ip1=i+1;
    
    #ifdef GLOBAL_MODEL
    if (i==nx-1){
		xP = (x-fx)/(2.*M_PI+ff[0].s0-fx); 
		ip1=0;
		if (!((xP<=1) && (xP>=0) && (xP==xP))){// in case of last phi=360 and first phi=0
			xP = (x-fx)/(2.*M_PI+ff[1].s0-fx); 
			ip1=1;
		}	
	}
	else{
	#endif
		xP = (x-fx)/(ff[i+1].s0-fx);
	#ifdef GLOBAL_MODEL
		}
	#endif
    double yP = (y-fy)/(ff[j+1].s1-fy);
    double zP = (z-fz)/(ff[k+1].s2-fz);
    
    double x0 = 1.0 - xP;
    double y0 = 1.0 - yP;
    double z0 = 1.0 - zP;
    
    
    
    
    size_t j0 =nx*j;
    size_t jp1=nx*(j+1);

    size_t k0 =nx*ny*k;
    size_t kp1=nx*ny*(k+1);

    
    double4 ind000 = B[i  +j0 +k0 ];
    double4 indP00 = B[ip1+j0 +k0 ];
    double4 ind0P0 = B[i  +jp1+k0 ];
    double4 ind00P = B[i  +j0 +kp1];
    double4 ind0PP = B[i  +jp1+kp1];
    double4 indP0P = B[ip1+j0 +kp1];
    double4 indPP0 = B[ip1+jp1+k0 ];
    double4 indPPP = B[ip1+jp1+kp1];
    
    double4 aa=a*( ind000*x0*y0*z0 + 
                   indP00*xP*y0*z0 + 
                   ind0P0*x0*yP*z0 + 
                   ind00P*x0*y0*zP + 
                   ind0PP*x0*yP*zP + 
                   indP0P*xP*y0*zP + 
                   indPP0*xP*yP*z0 + 
                   indPPP*xP*yP*zP) ;
	if (ind.s3 != 0) 
		aa.s3=sqrt(aa.s0*aa.s0+aa.s1*aa.s1+aa.s2*aa.s2);
	else
		aa.s3=0;
    return aa;
    

);

#if CALCULATE!=QSL
VEX_FUNCTION(cl_double4, interp_trilinearMid, (double, x)(double, y)(double,z)(cl_double4*, B)(cl_double4*, ff)(cl_ushort4, ind),
    DEFS();
    if (ind.s3 == 0) return 0.0;
    
    #ifdef PERIODIC_XY
		x=fmod(x-ff[0].s0+ff[nx-2].s0-ff[0].s0,ff[nx-2].s0-ff[0].s0)+ff[0].s0;
		y=fmod(y-ff[0].s1+ff[ny-2].s1-ff[0].s1,ff[ny-2].s1-ff[0].s1)+ff[0].s1;
	#else 
		#ifndef GLOBAL_MODEL
			if (x<=ff[0].s0) 
				x=ff[0].s0+1.e-4;
			
			if (x>=ff[nx-2].s0)  
				x=ff[nx-2].s0-1.e-4; 
		#else
			x=fmod((x+M_PI*2.0000001-ff[0].s0),(2.*M_PI))+ff[0].s0;
		#endif
		if (y<=ff[0].s1) 
			y=ff[0].s1+1.e-4; 
		if (y>=ff[ny-2].s1) 
			y=ff[ny-2].s1-1.e-4; 
	#endif
	
    if (z<=z_minimum) 
        z=z_minimum+1.e-4; 

    if (z<=ff[0].s2) 
        z=ff[0].s2+1.e-4; 
    
    if (z>=ff[nz-2].s2)  
        z=ff[nz-2].s2-1.e-4; 
    
    
    size_t i = ind.s0;
    size_t j = ind.s1;
    size_t k = ind.s2;
    
    double fx=ff[i].s0;
    double fy=ff[j].s1;
    double fz=ff[k].s2;
    
    double xP;
    size_t ip1=i+1;
    #ifdef GLOBAL_MODEL
    if (i==nx-2){
			xP = (x-fx)/(2.*M_PI+ff[0].s0-fx);
			ip1=0;	
		}
	else{
	#endif
		xP = (x-fx)/(ff[i+1].s0-fx);
	#ifdef GLOBAL_MODEL
		}
	#endif
    double yP = (y-fy)/(ff[j+1].s1-fy);
    double zP = (z-fz)/(ff[k+1].s2-fz);
    
    double x0 = 1.0 - xP;
    double y0 = 1.0 - yP;
    double z0 = 1.0 - zP;
    
    size_t j0 =nx*j;
    size_t jp1=nx*(j+1);

    size_t k0 =nx*ny*k;
    size_t kp1=nx*ny*(k+1);

    
    double4 ind000 = B[i  +j0 +k0 ];
    double4 indP00 = B[ip1+j0 +k0 ];
    double4 ind0P0 = B[i  +jp1+k0 ];
    double4 ind00P = B[i  +j0 +kp1];
    double4 ind0PP = B[i  +jp1+kp1];
    double4 indP0P = B[ip1+j0 +kp1];
    double4 indPP0 = B[ip1+jp1+k0 ];
    double4 indPPP = B[ip1+jp1+kp1];
    
    double4 aa= ind000*x0*y0*z0 + 
                   indP00*xP*y0*z0 + 
                   ind0P0*x0*yP*z0 + 
                   ind00P*x0*y0*zP + 
                   ind0PP*x0*yP*zP + 
                   indP0P*xP*y0*zP + 
                   indPP0*xP*yP*z0 + 
                   indPPP*xP*yP*zP ;
    
    return aa;
    

);
#endif


VEX_FUNCTION(cl_double4, interp_trilinear_diff, (double, x)(double, y)(double,z)(double, dx)(double, dy)\
             (double, dz)(cl_double4*, B)(cl_double4*, ff)(cl_ushort4, ind),
    if (ind.s3==0) return 0.0;
    DEFS();
    
    size_t i =ind.s0;
    size_t j =ind.s1;
    size_t k =ind.s2;
    
    double fx=ff[i].s0;
    double fy=ff[j].s1;
    double fz=ff[k].s2;

	#ifdef PERIODIC_XY
		x=fmod(x-ff[0].s0+ff[nx-1].s0-ff[0].s0,ff[nx-1].s0-ff[0].s0)+ff[0].s0;
		y=fmod(y-ff[0].s1+ff[ny-1].s1-ff[0].s1,ff[ny-1].s1-ff[0].s1)+ff[0].s1;
	#endif
	#ifdef GLOBAL_MODEL
		x=fmod((x+M_PI*2.0000001-ff[0].s0),(2.*M_PI))+ff[0].s0;
    #endif    

    double xP = (x-fx);
    double yP = (y-fy);
    double zP = (z-fz);
    
    
    size_t ip1=i+1;
    
    #ifdef GLOBAL_MODEL
    if (i==nx-1){
		fx-=2.*M_PI+ff[0].s0;
		ip1=0;
		xP /=-fx;
		if (!((xP<=1) && (xP>=0) && (xP==xP))){// in case of last phi=360 and first phi=0
			fx=-2.*M_PI-ff[1].s0+ff[i].s0;
			xP = (x-ff[i].s0)/(-fx); 
			ip1=1;
		}	
	}
	else{
	#endif
		fx-=ff[i+1].s0;
		xP /=-fx;
    #ifdef GLOBAL_MODEL
	}
    #endif
    fy-=ff[j+1].s1;
    fz-=ff[k+1].s2;
    
    
    yP /=-fy;
    zP /=-fz;
    
    double x0 = 1.0 - xP;
    double y0 = 1.0 - yP;
    double z0 = 1.0 - zP;
    
    double dxP = -dx/fx;
    double dyP = -dy/fy;
    double dzP = -dz/fz;
    
    double dx0 = -dxP;
    double dy0 = -dyP;
    double dz0 = -dzP;

    
    
    size_t j0 =nx*j;
    size_t jp1=nx*(j+1);

    size_t k0 =nx*ny*k;
    size_t kp1=nx*ny*(k+1);
    
    double4 ind000 = B[i  +j0 +k0 ];
    double4 indP00 = B[ip1+j0 +k0 ];
    double4 ind0P0 = B[i  +jp1+k0 ];
    double4 ind00P = B[i  +j0 +kp1];
    double4 ind0PP = B[i  +jp1+kp1];
    double4 indP0P = B[ip1+j0 +kp1];
    double4 indPP0 = B[ip1+jp1+k0 ];
    double4 indPPP = B[ip1+jp1+kp1];
    
    double4 aa = ind000*(dx0*y0*z0 + x0*dy0*z0 + x0*y0*dz0) + 
                 indP00*(dxP*y0*z0 + xP*dy0*z0 + xP*y0*dz0) + 
                 ind0P0*(dx0*yP*z0 + x0*dyP*z0 + x0*yP*dz0) + 
                 ind00P*(dx0*y0*zP + x0*dy0*zP + x0*y0*dzP) + 
                 ind0PP*(dx0*yP*zP + x0*dyP*zP + x0*yP*dzP) + 
                 indP0P*(dxP*y0*zP + xP*dy0*zP + xP*y0*dzP) + 
                 indPP0*(dxP*yP*z0 + xP*dyP*z0 + xP*yP*dz0) + 
                 indPPP*(dxP*yP*zP + xP*dyP*zP + xP*yP*dzP) ;

    return aa;
);
