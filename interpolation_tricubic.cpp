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

//// Used R. Keys, (1981). "Cubic convolution interpolation for digital image processing". IEEE Transactions on Acoustics, Speech, and Signal Processing 29 (6): 1153-1160. doi:10.1109/TASSP.1981.1163711.
//// No offset needed, contrary to the quad kernel, so x is positive.
#define W_cubic(x) \
    if (x<1.0) \
        x=1.5*x*x*x-2.5*x*x+1.0; \
    else if (x<2.0) \
        x=-0.5*x*x*x+2.5*x*x-4.0*x+2.0;  \
    else \
        x=0.0;

#define W_cubic_D(x) \
    if (x<1.0) \
        x=4.5*x*x-5.0*x; \
    else if (x<2.0) \
        x=-1.5*x*x+5.0*x-4.0;  \
    else \
        x=0.0;

////
VEX_FUNCTION(cl_double4, interp_tricubic,(double, x)(double, y)(double,z)(cl_double4*, B)(cl_double4*, ff)(cl_ushort4, ind),
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
    
    double fx=ff[i].s0;
    double fy=ff[j].s1;
    double fz=ff[k].s2;
    
    double x0 = (x-fx)/(ff[i+1].s0-fx);
    double y0 = (y-fy)/(ff[j+1].s1-fy);
    double z0 = (z-fz)/(ff[k+1].s2-fz);
    
    ////
    size_t im1=i-1;
    size_t ip1=i+1;
    size_t ip2=i+2;
    
    size_t j0 =nx*j;
    size_t jm1=nx*(j-1);
    size_t jp1=nx*(j+1);
    size_t jp2=nx*(j+2);
    
    size_t k0 =nx*ny*k;
    size_t km1=nx*ny*(k-1);
    size_t kp1=nx*ny*(k+1);
    size_t kp2=nx*ny*(k+2);
    
    
    
    
    double4 ind000 = B[i  +j0 +k0 ];
    double4 indP00 = B[ip1+j0 +k0 ];
    double4 ind0P0 = B[i  +jp1+k0 ];
    double4 ind00P = B[i  +j0 +kp1];
    double4 ind0PP = B[i  +jp1+kp1];
    double4 indP0P = B[ip1+j0 +kp1];
    double4 indPP0 = B[ip1+jp1+k0 ];
    double4 indPPP = B[ip1+jp1+kp1];
    
    double4 q;
    
    if ((i==0) || (j==0) || (k==0) || (i>nx-3) || (j>ny-3) || (k>nz-3)){
		    //Default to trilinear interpolation near edges.
            double xP=x0;
            double yP=y0;
            double zP=z0;
            
            x0 = 1.0 - xP;
            y0 = 1.0 - yP;
            z0 = 1.0 - zP;   
            q=a*( ind000*x0*y0*z0 + 
                   indP00*xP*y0*z0 + 
                   ind0P0*x0*yP*z0 + 
                   ind00P*x0*y0*zP + 
                   ind0PP*x0*yP*zP + 
                   indP0P*xP*y0*zP + 
                   indPP0*xP*yP*z0 + 
                   indPPP*xP*yP*zP) ;

        q.s3=sqrt(q.s0*q.s0+q.s1*q.s1+q.s2*q.s2);

            return q;
    }
    

    
    double x2 = 2.0 - x0;
    double y2 = 2.0 - y0;
    double z2 = 2.0 - z0;
    
    double xP = 1.0 - x0;
    double yP = 1.0 - y0;
    double zP = 1.0 - z0;
    
    double xM = 1.0 + x0;
    double yM = 1.0 + y0;
    double zM = 1.0 + z0;
 
 //////
    W_cubic(x0);
    W_cubic(y0);
    W_cubic(z0);
    
    W_cubic(xM);
    W_cubic(yM);
    W_cubic(zM);
    
    W_cubic(xP);
    W_cubic(yP);
    W_cubic(zP);
    
    W_cubic(x2);
    W_cubic(y2);
    W_cubic(z2);
    
    
    
//////
    double4 indM00 = B[im1+j0 +k0 ];
    double4 ind0M0 = B[i  +jm1+k0 ];
    double4 ind00M = B[i  +j0 +km1];
    double4 ind0MM = B[i  +jm1+km1];
    double4 indM0M = B[im1+j0 +km1];
    double4 indMM0 = B[im1+jm1+k0 ];
    double4 indMMM = B[im1+jm1+km1];

    double4 indMPP = B[im1+jp1+kp1];
    double4 indPMP = B[ip1+jm1+kp1];
    double4 indPPM = B[ip1+jp1+km1];
    double4 indPMM = B[ip1+jm1+km1];
    double4 indMPM = B[im1+jp1+km1];
    double4 indMMP = B[im1+jm1+kp1];

    double4 indMP0 = B[im1+jp1+k0 ];
    double4 indPM0 = B[ip1+jm1+k0 ];
    double4 indM0P = B[im1+j0 +kp1];
    double4 indP0M = B[ip1+j0 +km1];
    double4 ind0MP = B[i  +jm1+kp1];
    double4 ind0PM = B[i  +jp1+km1];

    double4 ind222 = B[ip2+jp2+kp2];

    double4 indP22 = B[ip1+jp2+kp2];
    double4 ind2P2 = B[ip2+jp1+kp2];
    double4 ind22P = B[ip2+jp2+kp1];
    double4 ind2PP = B[ip2+jp1+kp1];
    double4 indP2P = B[ip1+jp2+kp1];
    double4 indPP2 = B[ip1+jp1+kp2];

    double4 indM22 = B[im1+jp2+kp2];
    double4 ind2M2 = B[ip2+jm1+kp2];
    double4 ind22M = B[ip2+jp2+km1];
    double4 ind2MM = B[ip2+jm1+km1];
    double4 indM2M = B[im1+jp2+km1];
    double4 indMM2 = B[im1+jm1+kp2];

    double4 indMP2 = B[im1+jp1+kp2];
    double4 indPM2 = B[ip1+jm1+kp2];
    double4 indM2P = B[im1+jp2+kp1];
    double4 indP2M = B[ip1+jp2+km1];
    double4 ind2MP = B[ip2+jm1+kp1];
    double4 ind2PM = B[ip2+jp1+km1];

    double4 ind200 = B[ip2+j0 +k0 ];
    double4 ind020 = B[i  +jp2+k0 ];
    double4 ind002 = B[i  +j0 +kp2];
    double4 ind022 = B[i  +jp2+kp2];
    double4 ind202 = B[ip2+j0 +kp2];
    double4 ind220 = B[ip2+jp2+k0 ];

    double4 ind2P0 = B[ip2+jp1+k0 ];
    double4 indP20 = B[ip1+jp2+k0 ];
    double4 ind20P = B[ip2+j0 +kp1];
    double4 indP02 = B[ip1+j0 +kp2];
    double4 ind02P = B[i  +jp2+kp1];
    double4 ind0P2 = B[i  +jp1+kp2];

    double4 indM20 = B[im1+jp2+k0 ];
    double4 ind2M0 = B[ip2+jm1+k0 ];
    double4 indM02 = B[im1+j0 +kp2];
    double4 ind20M = B[ip2+j0 +km1];
    double4 ind0M2 = B[i  +jm1+kp2];
    double4 ind02M = B[i  +jp2+km1];
    
    q   =   a*(ind000*x0*y0*z0 + 
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
               ind0PM*x0*yP*zM +
               
               ind222*x2*y2*z2 +
               indP22*xP*y2*z2 +
               ind2P2*x2*yP*z2 +
               ind22P*x2*y2*zP +
               ind2PP*x2*yP*zP +
               indP2P*xP*y2*zP +
               indPP2*xP*yP*z2 +
               
               indM22*xM*y2*z2 +
               ind2M2*x2*yM*z2 +
               ind22M*x2*y2*zM +
               ind2MM*x2*yM*zM +
               indM2M*xM*y2*zM +
               indMM2*xM*yM*z2 +
                       
               indMP2*xM*yP*z2 +
               indPM2*xP*yM*z2 +
               indM2P*xM*y2*zP +
               indP2M*xP*y2*zM +
               ind2MP*x2*yM*zP +
               ind2PM*x2*yP*zM +
                       
               ind200*x2*y0*z0 +
               ind020*x0*y2*z0 +
               ind002*x0*y0*z2 +
               ind022*x0*y2*z2 +
               ind202*x2*y0*z2 +
               ind220*x2*y2*z0 +
                       
               ind2P0*x2*yP*z0 +
               indP20*xP*y2*z0 +
               ind20P*x2*y0*zP +
               indP02*xP*y0*z2 +
               ind02P*x0*y2*zP +
               ind0P2*x0*yP*z2 +
                       
               indM20*xM*y2*z0 +
               ind2M0*x2*yM*z0 +
               indM02*xM*y0*z2 +
               ind20M*x2*y0*zM +
               ind0M2*x0*yM*z2 +
               ind02M*x0*y2*zM );
           
        q.s3=sqrt(q.s0*q.s0+q.s1*q.s1+q.s2*q.s2);

    return q;
);




VEX_FUNCTION(cl_double4, interp_tricubic_diff,(double, x)(double, y)(double,z)(double, dx)(double, dy)\
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
    size_t ip2=i+2;
    
    size_t j0 =nx*j;
    size_t jm1=nx*(j-1);
    size_t jp1=nx*(j+1);
    size_t jp2=nx*(j+2);
    
    size_t k0 =nx*ny*k;
    size_t km1=nx*ny*(k-1);
    size_t kp1=nx*ny*(k+1);
    size_t kp2=nx*ny*(k+2);
    
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
    
            q =  ind000*(dx0*y0*z0 + x0*dy0*z0 + x0*y0*dz0) + 
                 indP00*(dxP*y0*z0 + xP*dy0*z0 + xP*y0*dz0) + 
                 ind0P0*(dx0*yP*z0 + x0*dyP*z0 + x0*yP*dz0) + 
                 ind00P*(dx0*y0*zP + x0*dy0*zP + x0*y0*dzP) + 
                 ind0PP*(dx0*yP*zP + x0*dyP*zP + x0*yP*dzP) + 
                 indP0P*(dxP*y0*zP + xP*dy0*zP + xP*y0*dzP) + 
                 indPP0*(dxP*yP*z0 + xP*dyP*z0 + xP*yP*dz0) + 
                 indPPP*(dxP*yP*zP + xP*dyP*zP + xP*yP*dzP) ;

            return q;
    }
    
    double x2 = 2.0 - x0;
    double y2 = 2.0 - y0;
    double z2 = 2.0 - z0;
    
    double xP = 1.0 - x0;
    double yP = 1.0 - y0;
    double zP = 1.0 - z0;
    
    double xM = 1.0 + x0;
    double yM = 1.0 + y0;
    double zM = 1.0 + z0;
    
    double dx0 = x0;
    double dy0 = y0;
    double dz0 = z0;
                   
    double dxP = xP;
    double dyP = yP;
    double dzP = zP;
                   
                   
    double dxM = xM;
    double dyM = yM;
    double dzM = zM;
    
    double dx2 = x2;
    double dy2 = y2;
    double dz2 = z2;
 //////
    W_cubic(x0);
    W_cubic(y0);
    W_cubic(z0);
    
    W_cubic(xM);
    W_cubic(yM);
    W_cubic(zM);
    
    W_cubic(xP);
    W_cubic(yP);
    W_cubic(zP);
    
    W_cubic(x2);
    W_cubic(y2);
    W_cubic(z2);
    
 //////
    W_cubic_D(dx0);
    W_cubic_D(dy0);
    W_cubic_D(dz0);
        
    W_cubic_D(dxP);
    W_cubic_D(dyP);
    W_cubic_D(dzP);
        
    W_cubic_D(dxM);
    W_cubic_D(dyM);
    W_cubic_D(dzM);
      
    W_cubic_D(dx2);
    W_cubic_D(dy2);
    W_cubic_D(dz2);
    
    dx0 *= dx*fx;
    dy0 *= dy*fy;
    dz0 *= dz*fz;
        
    dxP *= -dx*fx;
    dyP *= -dy*fy;
    dzP *= -dz*fz;
        
        
    dxM *= dx*fx;
    dyM *= dy*fy;
    dzM *= dz*fz;
    
    dx2 *= -dx*fx;
    dy2 *= -dy*fy;
    dz2 *= -dz*fz;
    
    
    
////
//////
    double4 indM00 = B[im1+j0 +k0 ];
    double4 ind0M0 = B[i  +jm1+k0 ];
    double4 ind00M = B[i  +j0 +km1];
    double4 ind0MM = B[i  +jm1+km1];
    double4 indM0M = B[im1+j0 +km1];
    double4 indMM0 = B[im1+jm1+k0 ];
    double4 indMMM = B[im1+jm1+km1];

    double4 indMPP = B[im1+jp1+kp1];
    double4 indPMP = B[ip1+jm1+kp1];
    double4 indPPM = B[ip1+jp1+km1];
    double4 indPMM = B[ip1+jm1+km1];
    double4 indMPM = B[im1+jp1+km1];
    double4 indMMP = B[im1+jm1+kp1];

    double4 indMP0 = B[im1+jp1+k0 ];
    double4 indPM0 = B[ip1+jm1+k0 ];
    double4 indM0P = B[im1+j0 +kp1];
    double4 indP0M = B[ip1+j0 +km1];
    double4 ind0MP = B[i  +jm1+kp1];
    double4 ind0PM = B[i  +jp1+km1];

    double4 ind222 = B[ip2+jp2+kp2];

    double4 indP22 = B[ip1+jp2+kp2];
    double4 ind2P2 = B[ip2+jp1+kp2];
    double4 ind22P = B[ip2+jp2+kp1];
    double4 ind2PP = B[ip2+jp1+kp1];
    double4 indP2P = B[ip1+jp2+kp1];
    double4 indPP2 = B[ip1+jp1+kp2];

    double4 indM22 = B[im1+jp2+kp2];
    double4 ind2M2 = B[ip2+jm1+kp2];
    double4 ind22M = B[ip2+jp2+km1];
    double4 ind2MM = B[ip2+jm1+km1];
    double4 indM2M = B[im1+jp2+km1];
    double4 indMM2 = B[im1+jm1+kp2];

    double4 indMP2 = B[im1+jp1+kp2];
    double4 indPM2 = B[ip1+jm1+kp2];
    double4 indM2P = B[im1+jp2+kp1];
    double4 indP2M = B[ip1+jp2+km1];
    double4 ind2MP = B[ip2+jm1+kp1];
    double4 ind2PM = B[ip2+jp1+km1];

    double4 ind200 = B[ip2+j0 +k0 ];
    double4 ind020 = B[i  +jp2+k0 ];
    double4 ind002 = B[i  +j0 +kp2];
    double4 ind022 = B[i  +jp2+kp2];
    double4 ind202 = B[ip2+j0 +kp2];
    double4 ind220 = B[ip2+jp2+k0 ];

    double4 ind2P0 = B[ip2+jp1+k0 ];
    double4 indP20 = B[ip1+jp2+k0 ];
    double4 ind20P = B[ip2+j0 +kp1];
    double4 indP02 = B[ip1+j0 +kp2];
    double4 ind02P = B[i  +jp2+kp1];
    double4 ind0P2 = B[i  +jp1+kp2];

    double4 indM20 = B[im1+jp2+k0 ];
    double4 ind2M0 = B[ip2+jm1+k0 ];
    double4 indM02 = B[im1+j0 +kp2];
    double4 ind20M = B[ip2+j0 +km1];
    double4 ind0M2 = B[i  +jm1+kp2];
    double4 ind02M = B[i  +jp2+km1];
    
    q =       (ind000*(dx0*y0*z0   +   x0*dy0*z0   +   x0*y0*dz0 )+ 
               indP00*(dxP*y0*z0   +   xP*dy0*z0   +   xP*y0*dz0 )+ 
               ind0P0*(dx0*yP*z0   +   x0*dyP*z0   +   x0*yP*dz0 )+ 
               ind00P*(dx0*y0*zP   +   x0*dy0*zP   +   x0*y0*dzP )+ 
               ind0PP*(dx0*yP*zP   +   x0*dyP*zP   +   x0*yP*dzP )+ 
               indP0P*(dxP*y0*zP   +   xP*dy0*zP   +   xP*y0*dzP )+ 
               indPP0*(dxP*yP*z0   +   xP*dyP*z0   +   xP*yP*dz0 )+ 
               indPPP*(dxP*yP*zP   +   xP*dyP*zP   +   xP*yP*dzP )+
                   
               indM00*(dxM*y0*z0   +   xM*dy0*z0   +   xM*y0*dz0 )+
               ind0M0*(dx0*yM*z0   +   x0*dyM*z0   +   x0*yM*dz0 )+
               ind00M*(dx0*y0*zM   +   x0*dy0*zM   +   x0*y0*dzM )+
               ind0MM*(dx0*yM*zM   +   x0*dyM*zM   +   x0*yM*dzM )+
               indM0M*(dxM*y0*zM   +   xM*dy0*zM   +   xM*y0*dzM )+
               indMM0*(dxM*yM*z0   +   xM*dyM*z0   +   xM*yM*dz0 )+
               indMMM*(dxM*yM*zM   +   xM*dyM*zM   +   xM*yM*dzM )+
     
               indMPP*(dxM*yP*zP   +   xM*dyP*zP   +   xM*yP*dzP )+
               indPMP*(dxP*yM*zP   +   xP*dyM*zP   +   xP*yM*dzP )+
               indPPM*(dxP*yP*zM   +   xP*dyP*zM   +   xP*yP*dzM )+
               indPMM*(dxP*yM*zM   +   xP*dyM*zM   +   xP*yM*dzM )+
               indMPM*(dxM*yP*zM   +   xM*dyP*zM   +   xM*yP*dzM )+
               indMMP*(dxM*yM*zP   +   xM*dyM*zP   +   xM*yM*dzP )+

               indMP0*(dxM*yP*z0   +   xM*dyP*z0   +   xM*yP*dz0 )+
               indPM0*(dxP*yM*z0   +   xP*dyM*z0   +   xP*yM*dz0 )+
               indM0P*(dxM*y0*zP   +   xM*dy0*zP   +   xM*y0*dzP )+
               indP0M*(dxP*y0*zM   +   xP*dy0*zM   +   xP*y0*dzM )+
               ind0MP*(dx0*yM*zP   +   x0*dyM*zP   +   x0*yM*dzP )+
               ind0PM*(dx0*yP*zM   +   x0*dyP*zM   +   x0*yP*dzM )+

               ind222*(dx2*y2*z2   +   x2*dy2*z2   +   x2*y2*dz2 )+
               indP22*(dxP*y2*z2   +   xP*dy2*z2   +   xP*y2*dz2 )+
               ind2P2*(dx2*yP*z2   +   x2*dyP*z2   +   x2*yP*dz2 )+
               ind22P*(dx2*y2*zP   +   x2*dy2*zP   +   x2*y2*dzP )+
               ind2PP*(dx2*yP*zP   +   x2*dyP*zP   +   x2*yP*dzP )+
               indP2P*(dxP*y2*zP   +   xP*dy2*zP   +   xP*y2*dzP )+
               indPP2*(dxP*yP*z2   +   xP*dyP*z2   +   xP*yP*dz2 )+
                   
               indM22*(dxM*y2*z2   +   xM*dy2*z2   +   xM*y2*dz2 )+
               ind2M2*(dx2*yM*z2   +   x2*dyM*z2   +   x2*yM*dz2 )+
               ind22M*(dx2*y2*zM   +   x2*dy2*zM   +   x2*y2*dzM )+
               ind2MM*(dx2*yM*zM   +   x2*dyM*zM   +   x2*yM*dzM )+
               indM2M*(dxM*y2*zM   +   xM*dy2*zM   +   xM*y2*dzM )+
               indMM2*(dxM*yM*z2   +   xM*dyM*z2   +   xM*yM*dz2 )+

               indMP2*(dxM*yP*z2   +   xM*dyP*z2   +   xM*yP*dz2 )+
               indPM2*(dxP*yM*z2   +   xP*dyM*z2   +   xP*yM*dz2 )+
               indM2P*(dxM*y2*zP   +   xM*dy2*zP   +   xM*y2*dzP )+
               indP2M*(dxP*y2*zM   +   xP*dy2*zM   +   xP*y2*dzM )+
               ind2MP*(dx2*yM*zP   +   x2*dyM*zP   +   x2*yM*dzP )+
               ind2PM*(dx2*yP*zM   +   x2*dyP*zM   +   x2*yP*dzM )+

               ind200*(dx2*y0*z0   +   x2*dy0*z0   +   x2*y0*dz0 )+
               ind020*(dx0*y2*z0   +   x0*dy2*z0   +   x0*y2*dz0 )+
               ind002*(dx0*y0*z2   +   x0*dy0*z2   +   x0*y0*dz2 )+
               ind022*(dx0*y2*z2   +   x0*dy2*z2   +   x0*y2*dz2 )+
               ind202*(dx2*y0*z2   +   x2*dy0*z2   +   x2*y0*dz2 )+
               ind220*(dx2*y2*z0   +   x2*dy2*z0   +   x2*y2*dz0 )+

               ind2P0*(dx2*yP*z0   +   x2*dyP*z0   +   x2*yP*dz0 )+
               indP20*(dxP*y2*z0   +   xP*dy2*z0   +   xP*y2*dz0 )+
               ind20P*(dx2*y0*zP   +   x2*dy0*zP   +   x2*y0*dzP )+
               indP02*(dxP*y0*z2   +   xP*dy0*z2   +   xP*y0*dz2 )+
               ind02P*(dx0*y2*zP   +   x0*dy2*zP   +   x0*y2*dzP )+
               ind0P2*(dx0*yP*z2   +   x0*dyP*z2   +   x0*yP*dz2 )+

               indM20*(dxM*y2*z0   +   xM*dy2*z0   +   xM*y2*dz0 )+
               ind2M0*(dx2*yM*z0   +   x2*dyM*z0   +   x2*yM*dz0 )+
               indM02*(dxM*y0*z2   +   xM*dy0*z2   +   xM*y0*dz2 )+
               ind20M*(dx2*y0*zM   +   x2*dy0*zM   +   x2*y0*dzM )+
               ind0M2*(dx0*yM*z2   +   x0*dyM*z2   +   x0*yM*dz2 )+
               ind02M*(dx0*y2*zM   +   x0*dy2*zM   +   x0*y2*dzM ));
                   
    return q;
            
);   
