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
void find_hft(void){
    using namespace std;
    double norm,trace;
    
    ofstream myfile0,myfile1,myfile2;
    ofstream myfile3,myfile4;
    #if RHO_Z==SYMM_LAMBDA
      myfile0.open ("ReDeltaLambda_symm.dat");
    #elif RHO_Z==OPT2_LAMBDA
      myfile0.open ("ReDeltaLambda_option2.dat");
    #else
      myfile0.open ("ReDeltaLambda.dat");
    #endif

    myfile1.open ("ODE_type.dat");
    myfile2.open ("Alpha_im.dat");
    myfile3.open ("ImLambda.dat");
    myfile4.open ("Trace.dat");
    
    
    cerr << "Analyzing transverse FL motions ...\n";
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
            //use trilinear interpolation to find grad B at cell center.
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
            norm=pow(bx000*bx000+by000*by000+bz000*bz000,0.5)+1.e-30;
            bx000/=norm;
            by000/=norm;
            bz000/=norm;
            
            double bx001  = Bx[i+nx*(j+ny*(k+1))];
            double by001  = By[i+nx*(j+ny*(k+1))];
            double bz001  = Bz[i+nx*(j+ny*(k+1))];
            norm=pow(bx001*bx001+by001*by001+bz001*bz001,0.5)+1.e-30;
            bx001/=norm;
            by001/=norm;
            bz001/=norm;
             
            double bx010  = Bx[i+nx*(j+1+ny*(k))];
            double by010  = By[i+nx*(j+1+ny*(k))];
            double bz010  = Bz[i+nx*(j+1+ny*(k))];
            norm=pow(bx010*bx010+by010*by010+bz010*bz010,0.5)+1.e-30;
            bx010/=norm;
            by010/=norm;
            bz010/=norm;
            
            double bx100  = Bx[i+1+nx*(j+ny*(k))];
            double by100  = By[i+1+nx*(j+ny*(k))];
            double bz100  = Bz[i+1+nx*(j+ny*(k))];
            norm=pow(bx100*bx100+by100*by100+bz100*bz100,0.5)+1.e-30;
            bx100/=norm;
            by100/=norm;
            bz100/=norm;
            
            double bx110  = Bx[i+1+nx*(j+1+ny*(k))];
            double by110  = By[i+1+nx*(j+1+ny*(k))];
            double bz110  = Bz[i+1+nx*(j+1+ny*(k))];
            norm=pow(bx110*bx110+by110*by110+bz110*bz110,0.5)+1.e-30;
            bx110/=norm;
            by110/=norm;
            bz110/=norm;
            
            double bx101  = Bx[i+1+nx*(j+ny*(k+1))];
            double by101  = By[i+1+nx*(j+ny*(k+1))];
            double bz101  = Bz[i+1+nx*(j+ny*(k+1))];
            norm=pow(bx101*bx101+by101*by101+bz101*bz101,0.5)+1.e-30;
            bx101/=norm;
            by101/=norm;
            bz101/=norm;
            
            double bx011  = Bx[i+nx*(j+1+ny*(k+1))];
            double by011  = By[i+nx*(j+1+ny*(k+1))];
            double bz011  = Bz[i+nx*(j+1+ny*(k+1))];
            norm=pow(bx011*bx011+by011*by011+bz011*bz011,0.5)+1.e-30;
            bx011/=norm;
            by011/=norm;
            bz011/=norm;
            
            double bx111  = Bx[i+1+nx*(j+1+ny*(k+1))];
            double by111  = By[i+1+nx*(j+1+ny*(k+1))];
            double bz111  = Bz[i+1+nx*(j+1+ny*(k+1))];
            norm=pow(bx111*bx111+by111*by111+bz111*bz111,0.5)+1.e-30;
            bx111/=norm;
            by111/=norm;
            bz111/=norm;
            
            
            double aX=bx000;
            double bX=(bx100-bx000)                                    /sx;
            double cX=(bx010-bx000)                                    /sy;
            double dX=(bx001-bx000)                                    /sz;
            double eX=(bx110-bx100-bx010+bx000)                        /sx/sy;
            double fX=(bx101-bx100-bx001+bx000)                        /sx/sz;
            double gX=(bx011-bx010-bx001+bx000)                        /sy/sz;
            double hX=(bx111-bx110-bx101-bx011+bx100+bx010+bx001-bx000)/sx/sy/sz;
             
            double ax =aX+(bX*sx+cX*sy+dX*sz)*0.5+0.25*(eX*sx*sy+fX*sx*sz+gX*sy*sz)+0.125*hX*sx*sy*sz;
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
              
            double ay =aY+(bY*sx+cY*sy+dY*sz)*0.5+0.25*(eY*sx*sy+fY*sx*sz+gY*sy*sz)+0.125*hY*sx*sy*sz;
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
              
            double az =aZ+(bZ*sx+cZ*sy+dZ*sz)*0.5+0.25*(eZ*sx*sy+fZ*sx*sz+gZ*sy*sz)+0.125*hZ*sx*sy*sz;
            double bz1=bZ+0.5*(eZ*sy+fZ*sz)+0.25*hZ*sy*sz;
            double bz2=cZ+0.5*(eZ*sx+gZ*sz)+0.25*hZ*sx*sz;
            double bz3=dZ+0.5*(gZ*sy+fZ*sx)+0.25*hZ*sx*sy;

            
            norm=pow(ax*ax+ay*ay+az*az,0.5)+1.e-30;
            ax/=norm;
            ay/=norm;
            az/=norm;
//************

		#if GEOMETRY==SPHERICAL
            double M[3][3]={	{	bx1+1./radM*az-tan(latM)/radM*ay, 	//B^{\hat \phi}_{\ ;\hat \phi}
									bx2,						  		//B^{\hat \phi}_{\ ;\hat \theta}
									bx3									//B^{\hat \phi}_{\ ;\hat r}
								},					  					
								{	by1+tan(latM)/radM*ax,			  	//B^{\hat \theta}_{\ ;\hat \phi}
									by2+1./radM*az,               		//B^{\hat \theta}_{\ ;\hat \theta}
									by3									//B^{\hat \theta}_{\ ;\hat r}
								},                    				 	
								{	bz1-1./radM*ax,					  	//B^{\hat r}_{\ ;\hat \phi}
									bz2-1./radM*ay,                		//B^{\hat r}_{\ ;\hat \theta}
									bz3									//B^{\hat r}_{\ ;\hat r}
									}
							};                    
       #else
            double M[3][3]={{bx1,bx2,bx3},
							{by1,by2,by3},
							{bz1,bz2,bz3}};
       #endif
        
            double alpha=IntegrandC[i+nx*(j+ny*k)];
            
            double P[3][3]={{1.0-ax*ax,-ax*ay,-ax*az},{-ax*ay,1.0-ay*ay,-ay*az},{-ax*az,-ay*az,1.0-az*az}};
            double J[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            for (uint r=0;r<3;r++)
                for (uint s=0;s<3;s++)
                    for (uint t=0;t<3;t++)
                        for (uint u=0;u<3;u++)
                            J[r][s]+=P[r][t]*M[t][u]*P[u][s]; //Take transverse piece of M
            
            trace=0.0;
            for (uint r=0;r<3;r++)
                        trace+=J[r][r];
            
            double traceM2=0;
            double traceMMt=0;
            for (uint r=0;r<3;r++)
                for (uint s=0;s<3;s++){
                        traceM2+=J[r][s]*J[s][r];
                        traceMMt+=pow(J[r][s],2);
					}
            
            double S=2.0*traceM2-trace*trace;
            //eigenvalues = (trace+/-sqrt(S))/2.0
            
            int Type;
            double dLambda,imLambda,dLambdaS;
            double curvature=0.;
            double curvfll=0.;
            
            dLambdaS=sqrt(traceMMt-(trace*trace-traceM2));//symmetrized lambda, always >= 0
            
            if (fabs(S)<trace*trace*1e-16){//repeated real roots
				int rank=2;
                if (fabs(traceMMt-traceM2)>trace*trace*1e-16) // works for 2x2 matrix
					rank=1; 
					
                if (rank==1)
                    Type=6;// improper node; one (repeated) eigenvector
                else
                    Type=7;// star node; two distinct eigenvectors
                
                dLambda=0.0;
                imLambda=0.0;
                IntegrandD[i+nx*(j+ny*k)]=0;
            }
            else if (S<0){//complex eigenvalues
                if (alpha>=0){
                    Type=3;//RH spiral
                    if (fabs(S)*1e-16>trace*trace)
                        Type=5;//RH center
                    imLambda=sqrt(fabs(S))/2.0;
                }
                else{
                    Type=2;//LH spiral
                    if (fabs(S)*1e-16>trace*trace)
                        Type=4;//LH center
                    imLambda=-sqrt(fabs(S))/2.0;
                }
                dLambda=0.0;
		#if RHO_Z==OPT2_LAMBDA
                	dLambda=2.*(2.*dLambdaS*dLambdaS-S)/(-S);
                	dLambda=log(dLambda/2.)/(M_PI_2/fabs(imLambda));
		#endif
            }
            else{
                double r1=(trace+sqrt(S))/2.0;
                double r2=(trace-sqrt(S))/2.0;
                IntegrandD[i+nx*(j+ny*k)]=0;
                if (r1*r2<=0){
                    Type=0; // saddle
				}
                else{
                    Type=1; // node
				}
				dLambda=sqrt(S);
				imLambda=0.0;
            }
		 //((m*m.T).trace()-(m.trace()^2-(m^2).trace()))
 		 #if RHO_Z==SYMM_LAMBDA 
		 	dLambda=dLambdaS
		 #endif
		 myfile0 << dLambda  << "\n";
		 myfile1 << Type 	 << "\n";
		 myfile2 << IntegrandD[i+nx*(j+ny*k)] 	 << "\n";
		 myfile3 << imLambda << "\n";
		 myfile4 << trace << "\n";
		 
		IntegrandA[i+nx*(j+ny*(k))]=imLambda;
		IntegrandB[i+nx*(j+ny*(k))]=dLambda;
		 
		 
       }
       myfile0.close();
       myfile1.close();
       myfile2.close();
       myfile3.close();
       myfile4.close();
       cerr << "... done.\n"; 
       
}
