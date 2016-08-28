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


////https://www.cl.cam.ac.uk/~nad10/pubs/quad.pdf
////NEEDS TO BE OFFSET BY 0.5. See interp_triquadratic for example. 
////Quadratic Interpolation for Image Resampling, Neil Dodgson, 1997, IEEE Transactions on Image Processing
//// r=1 to pass through points.
#define r_quad 1.0
#define W_quad(x) \
    x=fabs(x);\
    if (x<0.5) \
        x=-2.0*r_quad*x*x+0.5*(r_quad+1.0); \
    else if (x<1.5) \
        x=r_quad*x*x+(-2.0*r_quad-0.5)*x+0.75*(r_quad+1.0);  \
    else \
        x=0.0;
////
//// Note that due to offset of 0.5, x can be negative. So, take care of that below when taking derivative.
////
#define W_quad_D(x) \
    if (fabs(x)<0.5) \
        x=-4.0*r_quad*x; \
    else if (fabs(x)<1.5 && x>0) \
        x=2.0*r_quad*x+(-2.0*r_quad-0.5);  \
    else if (fabs(x)<1.5 && x<0) \
        x=2.0*r_quad*x-(-2.0*r_quad-0.5);  \
    else \
        x=0.0;
////
VEX_FUNCTION(cl_double4, interp_triquadratic, (double, x)(double, y)(double,z)(cl_double4*, B)(cl_double4*, ff)(cl_ushort4, ind),
        DEFS();
    double a=1.0;
    if (ind.s3 == 0) a=1e-9;
    


    if (x<ff[0].s0) 
        x=ff[0].s0+1.e-6;

    if (y<ff[0].s1) 
        y=ff[0].s1+1.e-6; 
        
    if (z<z_minimum) 
        z=z_minimum+1.e-6; 
    
    
    if (x>=ff[nx-1].s0)  
        x=ff[nx-1].s0-1.e-6; 
        
    if (y>=ff[ny-1].s1) 
        y=ff[ny-1].s1-1.e-6; 
        
    if (z>=ff[nz-1].s2)  
        z=ff[nz-1].s2-1.e-6; 

    
    
    
    size_t i = ind.s0;
    size_t j = ind.s1;
    size_t k = ind.s2;
    
    double x0 = (x-ff[i].s0)/(ff[i+1].s0-ff[i].s0);
    double y0 = (y-ff[j].s1)/(ff[j+1].s1-ff[j].s1);
    double z0 = (z-ff[k].s2)/(ff[k+1].s2-ff[k].s2);
    
    
    
    
    size_t im1=i-1;
    size_t ip1=i+1;
    
    size_t j0 =nx*j;
    size_t jm1=nx*(j-1);
    size_t jp1=nx*(j+1);

    size_t k0 =nx*ny*k;
    size_t km1=nx*ny*(k-1);
    size_t kp1=nx*ny*(k+1);
    
    double4 ind000 = B[i  +j0 +k0 ];
    double4 indP00 = B[ip1+j0 +k0 ];
    double4 ind0P0 = B[i  +jp1+k0 ];
    double4 ind00P = B[i  +j0 +kp1];
    double4 ind0PP = B[i  +jp1+kp1];
    double4 indP0P = B[ip1+j0 +kp1];
    double4 indPP0 = B[ip1+jp1+k0 ];
    double4 indPPP = B[ip1+jp1+kp1];
    
    double4 aa;
    
    if ((i==0) || (j==0) || (k==0)|| (i>nx-3) || (j>ny-3) || (k>nz-3)){
     		//Default to trilinear interpolation near edges.
            double xP=x0;
            double yP=y0;
            double zP=z0;
            
            x0 = 1.0 - xP;
            y0 = 1.0 - yP;
            z0 = 1.0 - zP;   
            aa=a*( ind000*x0*y0*z0 + 
                   indP00*xP*y0*z0 + 
                   ind0P0*x0*yP*z0 + 
                   ind00P*x0*y0*zP + 
                   ind0PP*x0*yP*zP + 
                   indP0P*xP*y0*zP + 
                   indPP0*xP*yP*z0 + 
                   indPPP*xP*yP*zP) ;

            aa.s3=sqrt(aa.s0*aa.s0+aa.s1*aa.s1+aa.s2*aa.s2);

            return aa;
    }
    
    double xP = 1.0 - x0;
    double yP = 1.0 - y0;
    double zP = 1.0 - z0;
    
    double xM = 1.0 + x0;
    double yM = 1.0 + y0;
    double zM = 1.0 + z0;
    
    if (x0>0.5){
		xM=x0;
		x0=x0-1.0;
		xP=1.0-x0;
		i=i+1;
		}
	if (y0>0.5){
		yM=y0;
		y0=y0-1.0;
		yP=1.0-y0;
		j=j+1;
		}
	if (z0>0.5){
		zM=z0;
		z0=z0-1.0;
		zP=1.0-z0;
		k=k+1;
		}    
		
    im1=i-1;
    ip1=i+1;
    
    j0 =nx*j;
    jm1=nx*(j-1);
    jp1=nx*(j+1);
    
    k0 =nx*ny*k;
    km1=nx*ny*(k-1);
    kp1=nx*ny*(k+1);
    
    ind000 = B[i  +j0 +k0 ];
    indP00 = B[ip1+j0 +k0 ];
    ind0P0 = B[i  +jp1+k0 ];
    ind00P = B[i  +j0 +kp1];
    ind0PP = B[i  +jp1+kp1];
    indP0P = B[ip1+j0 +kp1];
    indPP0 = B[ip1+jp1+k0 ];
    indPPP = B[ip1+jp1+kp1];
     
    W_quad(x0);
    W_quad(y0);
    W_quad(z0);
    
    W_quad(xP);
    W_quad(yP);
    W_quad(zP);
    
    W_quad(xM);
    W_quad(yM);
    W_quad(zM);

             
    double4 indM00 = B[im1+j0 +k0 ];
    double4 indMPP = B[im1+jp1+kp1];
    double4 indMP0 = B[im1+jp1+k0 ];
    double4 indM0P = B[im1+j0 +kp1];
     
    double4 ind0M0 = B[i  +jm1+k0 ];
    double4 indPMP = B[ip1+jm1+kp1];
    double4 indPM0 = B[ip1+jm1+k0 ];
    double4 ind0MP = B[i  +jm1+kp1];
     
    double4 ind00M = B[i  +j0 +km1];
    double4 indPPM = B[ip1+jp1+km1];
    double4 indP0M = B[ip1+j0 +km1];
    double4 ind0PM = B[i  +jp1+km1];
     
    double4 ind0MM = B[i  +jm1+km1];
    double4 indPMM = B[ip1+jm1+km1];
     
    double4 indM0M = B[im1+j0 +km1];
    double4 indMPM = B[im1+jp1+km1];
     
    double4 indMM0 = B[im1+jm1+k0 ];
    double4 indMMP = B[im1+jm1+kp1];
     
    double4 indMMM = B[im1+jm1+km1];
    

    aa = a*(ind000*x0*y0*z0 + 
            indP00*xP*y0*z0 + 
            ind0P0*x0*yP*z0 + 
            ind00P*x0*y0*zP + 
            ind0PP*x0*yP*zP + 
            indP0P*xP*y0*zP + 
            indPP0*xP*yP*z0 + 
            indPPP*xP*yP*zP +
                  
            indM00*xM*y0*z0 +
            ind0M0*x0*yM*z0 +
            ind00M*x0*y0*zM +
            ind0MM*x0*yM*zM +
            indM0M*xM*y0*zM +
            indMM0*xM*yM*z0 +
            indMMM*xM*yM*zM +
                  
            indMPP*xM*yP*zP +
            indPMP*xP*yM*zP +
            indPPM*xP*yP*zM +
            indPMM*xP*yM*zM +
            indMPM*xM*yP*zM +
            indMMP*xM*yM*zP +
                  
            indMP0*xM*yP*z0 +
            indPM0*xP*yM*z0 +
            indM0P*xM*y0*zP +
            indP0M*xP*y0*zM +
            ind0MP*x0*yM*zP +
            ind0PM*x0*yP*zM );

        aa.s3=sqrt(aa.s0*aa.s0+aa.s1*aa.s1+aa.s2*aa.s2);

    
    return aa;
            
);               
////
////
VEX_FUNCTION(cl_double4, interp_triquadratic_diff,(double, x)(double, y)(double,z)(double, dx)(double, dy)\
             (double, dz)(cl_double4*, B)(cl_double4*, ff)(cl_ushort4, ind),
    if (ind.s3==0) return 0.0;
    DEFS();
    
    
    
    size_t i =ind.s0;
    size_t j =ind.s1;
    size_t k =ind.s2;
    
    double fx = 1.0/(ff[i+1].s0-ff[i].s0);
    double fy = 1.0/(ff[j+1].s1-ff[j].s1);
    double fz = 1.0/(ff[k+1].s2-ff[k].s2);
    
    double x0 = (x-ff[i].s0)*fx;
    double y0 = (y-ff[j].s1)*fy;
    double z0 = (z-ff[k].s2)*fz;
    
    
    size_t im1=i-1;
    size_t ip1=i+1;
    
    size_t j0 =nx*j;
    size_t jm1=nx*(j-1);
    size_t jp1=nx*(j+1);

    size_t k0 =nx*ny*k;
    size_t km1=nx*ny*(k-1);
    size_t kp1=nx*ny*(k+1);

    
    double4 q;
    
    double4 ind000 = B[i  +j0 +k0 ];
    double4 indP00 = B[ip1+j0 +k0 ];
    double4 ind0P0 = B[i  +jp1+k0 ];
    double4 ind00P = B[i  +j0 +kp1];
    double4 ind0PP = B[i  +jp1+kp1];
    double4 indP0P = B[ip1+j0 +kp1];
    double4 indPP0 = B[ip1+jp1+k0 ];
    double4 indPPP = B[ip1+jp1+kp1];
          
    
    if ((i==0) || (j==0) || (k==0) || (i>nx-3) || (j>ny-3) || (k>nz-3)){
   	        //Default to trilinear interpolation near edges.
            double xP=x0;
            double yP=y0;
            double zP=z0;
            
            x0 = 1.0 - xP;
            y0 = 1.0 - yP;
            z0 = 1.0 - zP;   
            
            double dxP = dx*fx;
            double dyP = dy*fy;
            double dzP = dz*fz;
    
            double dx0 = -dxP;
            double dy0 = -dyP;
            double dz0 = -dzP;
    
            q = ind000*(dx0*y0*z0 + x0*dy0*z0 + x0*y0*dz0) + 
                indP00*(dxP*y0*z0 + xP*dy0*z0 + xP*y0*dz0) + 
                ind0P0*(dx0*yP*z0 + x0*dyP*z0 + x0*yP*dz0) + 
                ind00P*(dx0*y0*zP + x0*dy0*zP + x0*y0*dzP) + 
                ind0PP*(dx0*yP*zP + x0*dyP*zP + x0*yP*dzP) + 
                indP0P*(dxP*y0*zP + xP*dy0*zP + xP*y0*dzP) + 
                indPP0*(dxP*yP*z0 + xP*dyP*z0 + xP*yP*dz0) + 
                indPPP*(dxP*yP*zP + xP*dyP*zP + xP*yP*dzP) ;

            return q;
    }
    
    double xP = 1.0 - x0;
    double yP = 1.0 - y0;
    double zP = 1.0 - z0;
    
    double xM = 1.0 + x0;
    double yM = 1.0 + y0;
    double zM = 1.0 + z0;
    

 
    
    if (x0>0.5){
		xM=x0;
		x0=x0-1.0;
		xP=1.0-x0;
		i=i+1;
		}
	if (y0>0.5){
		yM=y0;
		y0=y0-1.0;
		yP=1.0-y0;
		j=j+1;
		}
	if (z0>0.5){
		zM=z0;
		z0=z0-1.0;
		zP=1.0-z0;
		k=k+1;
		}    
	double dx0 = x0;
    double dy0 = y0;
    double dz0 = z0;
                   
    double dxP = xP;
    double dyP = yP;
    double dzP = zP;
                   
                   
    double dxM = xM;
    double dyM = yM;
    double dzM = zM;
		
    im1=i-1;
    ip1=i+1;
    
    j0 =nx*j;
    jm1=nx*(j-1);
    jp1=nx*(j+1);
    
    k0 =nx*ny*k;
    km1=nx*ny*(k-1);
    kp1=nx*ny*(k+1);
    
    ind000 = B[i  +j0 +k0 ];
    indP00 = B[ip1+j0 +k0 ];
    ind0P0 = B[i  +jp1+k0 ];
    ind00P = B[i  +j0 +kp1];
    ind0PP = B[i  +jp1+kp1];
    indP0P = B[ip1+j0 +kp1];
    indPP0 = B[ip1+jp1+k0 ];
    indPPP = B[ip1+jp1+kp1];
     
 
 //////
    W_quad(x0);
    W_quad(y0);
    W_quad(z0);
    
    W_quad(xP);
    W_quad(yP);
    W_quad(zP);
    
    W_quad(xM);
    W_quad(yM);
    W_quad(zM);
 //////
    W_quad_D(dx0);
    W_quad_D(dy0);
    W_quad_D(dz0);
            
    W_quad_D(dxP);
    W_quad_D(dyP);
    W_quad_D(dzP);
            
    W_quad_D(dxM);
    W_quad_D(dyM);
    W_quad_D(dzM);
    
    dx0 *= dx*fx;
    dy0 *= dy*fy;
    dz0 *= dz*fz;
        
    dxP *= -dx*fx;
    dyP *= -dy*fy;
    dzP *= -dz*fz;
        
        
    dxM *= dx*fx;
    dyM *= dy*fy;
    dzM *= dz*fz;
    
    
    

    double4 indM00 = B[im1+j0 +k0 ];
    double4 indMPP = B[im1+jp1+kp1];
    double4 indMP0 = B[im1+jp1+k0 ];
    double4 indM0P = B[im1+j0 +kp1];

    double4 ind0M0 = B[i  +jm1+k0 ];
    double4 indPMP = B[ip1+jm1+kp1];
    double4 indPM0 = B[ip1+jm1+k0 ];
    double4 ind0MP = B[i  +jm1+kp1];

    double4 ind00M = B[i  +j0 +km1];
    double4 indPPM = B[ip1+jp1+km1];
    double4 indP0M = B[ip1+j0 +km1];
    double4 ind0PM = B[i  +jp1+km1];
 
    double4 ind0MM = B[i  +jm1+km1];
    double4 indPMM = B[ip1+jm1+km1];
     
    double4 indM0M = B[im1+j0 +km1];
    double4 indMPM = B[im1+jp1+km1];
     
    double4 indMM0 = B[im1+jm1+k0 ];
    double4 indMMP = B[im1+jm1+kp1];
     
    double4 indMMM = B[im1+jm1+km1];
    
                     
    
    
    q = (ind000*(dx0*y0*z0 +   x0*dy0*z0 + x0*y0*dz0) + 
         indP00*(dxP*y0*z0 +   xP*dy0*z0 + xP*y0*dz0) + 
         ind0P0*(dx0*yP*z0 +   x0*dyP*z0 + x0*yP*dz0) + 
         ind00P*(dx0*y0*zP +   x0*dy0*zP + x0*y0*dzP) + 
         ind0PP*(dx0*yP*zP +   x0*dyP*zP + x0*yP*dzP) + 
         indP0P*(dxP*y0*zP +   xP*dy0*zP + xP*y0*dzP) + 
         indPP0*(dxP*yP*z0 +   xP*dyP*z0 + xP*yP*dz0) + 
         indPPP*(dxP*yP*zP +   xP*dyP*zP + xP*yP*dzP) +
                                                    
         indM00*(dxM*y0*z0 +   xM*dy0*z0 + xM*y0*dz0) +
         ind0M0*(dx0*yM*z0 +   x0*dyM*z0 + x0*yM*dz0) +
         ind00M*(dx0*y0*zM +   x0*dy0*zM + x0*y0*dzM) +
         ind0MM*(dx0*yM*zM +   x0*dyM*zM + x0*yM*dzM) +
         indM0M*(dxM*y0*zM +   xM*dy0*zM + xM*y0*dzM) +
         indMM0*(dxM*yM*z0 +   xM*dyM*z0 + xM*yM*dz0) +
         indMMM*(dxM*yM*zM +   xM*dyM*zM + xM*yM*dzM) +
                
         indMPP*(dxM*yP*zP +   xM*dyP*zP + xM*yP*dzP) +
         indPMP*(dxP*yM*zP +   xP*dyM*zP + xP*yM*dzP) +
         indPPM*(dxP*yP*zM +   xP*dyP*zM + xP*yP*dzM) +
         indPMM*(dxP*yM*zM +   xP*dyM*zM + xP*yM*dzM) +
         indMPM*(dxM*yP*zM +   xM*dyP*zM + xM*yP*dzM) +
         indMMP*(dxM*yM*zP +   xM*dyM*zP + xM*yM*dzP) +
                
         indMP0*(dxM*yP*z0 +   xM*dyP*z0 + xM*yP*dz0) +
         indPM0*(dxP*yM*z0 +   xP*dyM*z0 + xP*yM*dz0) +
         indM0P*(dxM*y0*zP +   xM*dy0*zP + xM*y0*dzP) +
         indP0M*(dxP*y0*zM +   xP*dy0*zM + xP*y0*dzM) +
         ind0MP*(dx0*yM*zP +   x0*dyM*zP + x0*yM*dzP) +
         ind0PM*(dx0*yP*zM +   x0*dyP*zM + x0*yP*dzM)) ;
                   
    return q;
            
);   
