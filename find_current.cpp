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


#include <iostream>

void find_current(void){
    using namespace std;
    double norm,alpha,alphaHat;
    
    
    cerr << "Calculating current ...\n";
    ofstream myfile0,myfile1,myfile2,myfile3,myfile4;
    myfile0.open ("Jx.dat");
    myfile1.open ("Jy.dat");
    myfile2.open ("Jz.dat");
    myfile3.open ("Alpha.dat");
    
    
    for (uint64_t k=0;k<nz-1;k++)
		for (uint64_t j=0;j<ny-1;j++)
			for (uint64_t i=0;i<nx-1;i++)
        {

            double lon0=Hx[i];//if geometry is cartesian, lon,lat,rad correspond to x,y,z
            double lat0=Hy[j];
            double rad0=Hz[k];

            double lon1=Hx[i+1];
            double lat1=Hy[j+1];
            double rad1=Hz[k+1];
            
            double lonM=HxMid[i];
            double latM=HyMid[j];
            double radM=HzMid[k];
            
            //use trilinear interpolation to find J at cell center.
            
            #if GEOMETRY!=SPHERICAL
				double sx=-lon0;//-rad*lon*cos(lat);
				double sy=-lat0;//-rad*lat; // yl	
				double sz=-rad0;//-rad;
				sx+=lon1;//rad*lon*cos(lat);
				sy+=lat1;//rad*lat; // yl	
				sz+=rad1;//rad;
            #else
				double sx=-lon0*radM*cos(latM);//-rad*lon*cos(lat);
				double sy=-lat0*radM;//-rad*lat; // yl	
				double sz=-rad0;//-rad;
				sx+=lon1*radM*cos(latM);//rad*lon*cos(lat);
				sy+=lat1*radM;//rad*lat; // yl	
				sz+=rad1;//rad;
            #endif
            double bx000 = Bx[i+nx*(j+ny*(k))];
            double by000 = By[i+nx*(j+ny*(k))];
            double bz000 = Bz[i+nx*(j+ny*(k))];

            double bx001  = Bx[i+nx*(j+ny*(k+1))];
            double by001  = By[i+nx*(j+ny*(k+1))];
            double bz001  = Bz[i+nx*(j+ny*(k+1))];
            
            double bx010  = Bx[i+nx*(j+1+ny*(k))];
            double by010  = By[i+nx*(j+1+ny*(k))];
            double bz010  = Bz[i+nx*(j+1+ny*(k))];
           
            double bx100  = Bx[i+1+nx*(j+ny*(k))];
            double by100  = By[i+1+nx*(j+ny*(k))];
            double bz100  = Bz[i+1+nx*(j+ny*(k))];
            
            double bx110  = Bx[i+1+nx*(j+1+ny*(k))];
            double by110  = By[i+1+nx*(j+1+ny*(k))];
            double bz110  = Bz[i+1+nx*(j+1+ny*(k))];
            
            double bx101  = Bx[i+1+nx*(j+ny*(k+1))];
            double by101  = By[i+1+nx*(j+ny*(k+1))];
            double bz101  = Bz[i+1+nx*(j+ny*(k+1))];
            
            double bx011  = Bx[i+nx*(j+1+ny*(k+1))];
            double by011  = By[i+nx*(j+1+ny*(k+1))];
            double bz011  = Bz[i+nx*(j+1+ny*(k+1))];
            
            double bx111  = Bx[i+1+nx*(j+1+ny*(k+1))];
            double by111  = By[i+1+nx*(j+1+ny*(k+1))];
            double bz111  = Bz[i+1+nx*(j+1+ny*(k+1))];

            double aX=bx000;
            double bX=(bx100-bx000)                                    /sx;
            double cX=(bx010-bx000)                                    /sy;
            double dX=(bx001-bx000)                                    /sz;
            double eX=(bx110-bx100-bx010+bx000)                        /sx/sy;
            double fX=(bx101-bx100-bx001+bx000)                        /sx/sz;
            double gX=(bx011-bx010-bx001+bx000)                        /sy/sz;
            double hX=(bx111-bx110-bx101-bx011+bx100+bx010+bx001-bx000)/sx/sy/sz;
             
            double bx =aX+(bX*sx+cX*sy+dX*sz)*0.5+0.25*(eX*sx*sy+fX*sx*sz+gX*sy*sz)+0.125*hX*sx*sy*sz;
            double bx1=bX+0.5*(eX*sy+fX*sz)+0.25*hX*sy*sz;
            double bx2=cX+0.5*(eX*sx+gX*sz)+0.25*hX*sx*sz;
            double bx3=dX+0.5*(gX*sy+fX*sx)+0.25*hX*sx*sy;
             
            double aY=by000;                                           
            double bY=(by100-by000                                     )/sx;
            double cY=(by010-by000                                     )/sy;
            double dY=(by001-by000                                     )/sz;
            double eY=(by110-by100-by010+by000                         )/sx/sy;
            double fY=(by101-by100-by001+by000                         )/sx/sz;
            double gY=(by011-by010-by001+by000                         )/sy/sz;
            double hY=(by111-by110-by101-by011+by100+by010+by001-by000 )/sx/sy/sz;
              
            double by =aY+(bY*sx+cY*sy+dY*sz)*0.5+0.25*(eY*sx*sy+fY*sx*sz+gY*sy*sz)+0.125*hY*sx*sy*sz;
            double by1=bY+0.5*(eY*sy+fY*sz)+0.25*hY*sy*sz;
            double by2=cY+0.5*(eY*sx+gY*sz)+0.25*hY*sx*sz;
            double by3=dY+0.5*(gY*sy+fY*sx)+0.25*hY*sx*sy;
             
            double aZ=bz000;
            double bZ=(bz100-bz000                                     )/sx;
            double cZ=(bz010-bz000                                     )/sy;
            double dZ=(bz001-bz000                                     )/sz;
            double eZ=(bz110-bz100-bz010+bz000                         )/sx/sy;
            double fZ=(bz101-bz100-bz001+bz000                         )/sx/sz;
            double gZ=(bz011-bz010-bz001+bz000                         )/sy/sz;
            double hZ=(bz111-bz110-bz101-bz011+bz100+bz010+bz001-bz000 )/sx/sy/sz;
              
            double bz =aZ+(bZ*sx+cZ*sy+dZ*sz)*0.5+0.25*(eZ*sx*sy+fZ*sx*sz+gZ*sy*sz)+0.125*hZ*sx*sy*sz;
            double bz1=bZ+0.5*(eZ*sy+fZ*sz)+0.25*hZ*sy*sz;
            double bz2=cZ+0.5*(eZ*sx+gZ*sz)+0.25*hZ*sx*sz;
            double bz3=dZ+0.5*(gZ*sy+fZ*sx)+0.25*hZ*sx*sy;

//************

		#if GEOMETRY==SPHERICAL
            double M[3][3]={	{	bx1+1./radM*bz-tan(latM)/radM*by, 	//B^{\hat \phi}_{\ ;\hat \phi}
									bx2,						  		//B^{\hat \phi}_{\ ;\hat \theta}
									bx3									//B^{\hat \phi}_{\ ;\hat r}
								},					  					
								{	by1+tan(latM)/radM*bx,			  	//B^{\hat \theta}_{\ ;\hat \phi}
									by2+1./radM*bz,               		//B^{\hat \theta}_{\ ;\hat \theta}
									by3									//B^{\hat \theta}_{\ ;\hat r}
								},                    				 	
								{	bz1-1./radM*bx,					  	//B^{\hat r}_{\ ;\hat \phi}
									bz2-1./radM*by,                		//B^{\hat r}_{\ ;\hat \theta}
									bz3									//B^{\hat r}_{\ ;\hat r}
									}
							};                    
       #else
            double M[3][3]={{bx1,bx2,bx3},
							{by1,by2,by3},
							{bz1,bz2,bz3}};
       #endif
			
			double jx=M[2][1]-M[1][2];
			double jy=M[0][2]-M[2][0];
			double jz=M[1][0]-M[0][1];
			
				alpha=(jx*bx+jy*by+jz*bz)/(bx*bx+by*by+bz*bz+1e-60);
				myfile0 << jx << "\n";
				myfile1 << jy << "\n";
				myfile2 << jz << "\n";
				myfile3 << alpha << "\n";
				IntegrandC[i+nx*(j+ny*k)] = alpha;
				IntegrandD[i+nx*(j+ny*k)] = alpha; // will be modified in find_hft() to eliminate locations with real transverse eigenvalues
		}
		
       myfile0.close();
       myfile1.close();
       myfile2.close();
       myfile3.close();
       cerr << "... done.\n"; 
}
