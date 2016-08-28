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

void hilbert_to_coo(qsl_type &qsl,double *lon,double *lat,double *rad,size_t i){
    using namespace std;
	// Works both in cartesian and spherical coordinates
	// In cartesian coordinates, x=lon, y=lat, z=rad
    
	#if QSL_DIM==3
			uint64_t x,y,z;
            //if (!SAMPLER_INITIALIZED){
            //    ORIGIN_XY[1]/=solar_radius;
            //    ORIGIN_XY[0]/=solar_radius; 
            //    LX/=solar_radius;
            //    LY/=solar_radius; 
            //    SAMPLER_INITIALIZED=true;
            //}
            point( qsl[i].x, &x, &y, &z);
            double dx=(double(x))/(double(MAX_RES));
            double dy=(double(y))/(double(MAX_RES));
            double dz=(double(z))/(double(MAX_RES));
            
            *lon=dx*LX+ORIGIN_XY[0];
            *lat=dy*LY+ORIGIN_XY[1];
            z_sampler(dz);
            *rad=dz+solar_radius_sub;
	#endif
	#if QSL_DIM==2
     	    double xx,yy,zz;
			uint64_t x,y;
            if (!SAMPLER_INITIALIZED){
                // normalize vector (n) normal to slice
                double norm=pow(SLICE_NORMAL[0]*SLICE_NORMAL[0]+SLICE_NORMAL[1]*SLICE_NORMAL[1]+SLICE_NORMAL[2]*SLICE_NORMAL[2],0.5);
                SLICE_NORMAL[0]/=norm;
                SLICE_NORMAL[1]/=norm;
                SLICE_NORMAL[2]/=norm;
                // force u=SLICE_NORMAL to lie in the plane of the slice: u=u-n(u.n). Then normalize.
                double n_dot_u=SLICE_NORMAL[0]*SLICE_UP[0]+SLICE_NORMAL[1]*SLICE_UP[1]+SLICE_NORMAL[2]*SLICE_UP[2];
                SLICE_UP[0]-=SLICE_NORMAL[0]*n_dot_u;
                SLICE_UP[1]-=SLICE_NORMAL[1]*n_dot_u;
                SLICE_UP[2]-=SLICE_NORMAL[2]*n_dot_u;
                norm=pow(SLICE_UP[0]*SLICE_UP[0]+SLICE_UP[1]*SLICE_UP[1]+SLICE_UP[2]*SLICE_UP[2],0.5);
                SLICE_UP[0]/=norm;
                SLICE_UP[1]/=norm;
                SLICE_UP[2]/=norm;
                //do cross product between u and n to find second basis vector in the slice:
                SLICE_RIGHT[0]=SLICE_UP[1]*SLICE_NORMAL[2]-SLICE_UP[2]*SLICE_NORMAL[1];
                SLICE_RIGHT[1]=-SLICE_UP[0]*SLICE_NORMAL[2]+SLICE_UP[2]*SLICE_NORMAL[0];
                SLICE_RIGHT[2]=SLICE_UP[0]*SLICE_NORMAL[1]-SLICE_UP[1]*SLICE_NORMAL[0];
                //basis in plane of slice is constructed.

				SLICE_CENTER[0]*=to_radians;
				SLICE_CENTER[1]*=to_radians; 
                SLICE_CENTER[2]+=solar_radius_sub; // distance from center of sun for spherical coo. else this adds zero and does nothing.

                #if ((GEOMETRY==SPHERICAL) && (SLICE_TYPE==CARTESIAN))
					// Need to convert unit vectors from spherical coo to cartesian.
					// Use Jackson's conventions.
					double uf=  SLICE_UP[0];//phi
					double ut= -SLICE_UP[1];//theta
					double ur=  SLICE_UP[2];//r
					double rf=  SLICE_RIGHT[0];//phi
					double rt= -SLICE_RIGHT[1];//theta
					double rr=  SLICE_RIGHT[2];//r
					
					double ux,uy,uz,rx,ry,rz;
					double phi  = SLICE_CENTER[0];
					double theta= 3.1415926/2.0 - SLICE_CENTER[1];
					double r=SLICE_CENTER[2];
					
					ux=sin(theta)*cos(phi)*ur +
					   cos(theta)*cos(phi)*ut -
					   sin(phi)*uf;
					   
					uy=sin(theta)*sin(phi)*ur +
					   cos(theta)*sin(phi)*ut +
					   cos(phi)*uf;
					   
					uz=cos(theta)*ur -
					   sin(theta)*ut;
					
					rx=sin(theta)*cos(phi)*rr +
					   cos(theta)*cos(phi)*rt -
					   sin(phi)*rf;
					   
					ry=sin(theta)*sin(phi)*rr +
					   cos(theta)*sin(phi)*rt +
					   cos(phi)*rf;
					   
					rz=cos(theta)*rr -
					   sin(theta)*rt;
					
					SLICE_UP[0]  =ux; 
					SLICE_UP[1]  =uy; 
					SLICE_UP[2]  =uz;
					
					SLICE_RIGHT[0]  =rx; 
					SLICE_RIGHT[1]  =ry; 
					SLICE_RIGHT[2]  =rz;
					
					
					
					SLICE_CENTER[0]=r*sin(theta)*cos(phi);
					SLICE_CENTER[1]=r*sin(theta)*sin(phi);
					SLICE_CENTER[2]=r*cos(theta);
					
                #endif
                #if ((GEOMETRY==SPHERICAL) && (SLICE_TYPE==SPHERICAL))
					SLICE_LX*=to_radians; 
					SLICE_LY*=to_radians; 
				#endif
                
                //move SLICE_CENTER to be at origin of slice:
                
                SLICE_CENTER[0]-=(SLICE_LX/2.0)*SLICE_RIGHT[0]+(SLICE_LY/2.0)*SLICE_UP[0];
                SLICE_CENTER[1]-=(SLICE_LX/2.0)*SLICE_RIGHT[1]+(SLICE_LY/2.0)*SLICE_UP[1];
                SLICE_CENTER[2]-=(SLICE_LX/2.0)*SLICE_RIGHT[2]+(SLICE_LY/2.0)*SLICE_UP[2];
                
                
                
                SAMPLER_INITIALIZED=true;
            }
            
            
            point( qsl[i].x, &x, &y);
            
            double dx=((double)x)/((double)MAX_RES);
            double dy=((double)y)/((double)MAX_RES);
            
            xx=SLICE_CENTER[0];
            yy=SLICE_CENTER[1];
            zz=SLICE_CENTER[2];
            
            xx+=SLICE_LX*dx*SLICE_RIGHT[0];
            yy+=SLICE_LX*dx*SLICE_RIGHT[1];
            zz+=SLICE_LX*dx*SLICE_RIGHT[2];
            
            xx+=SLICE_LY*dy*SLICE_UP[0];
            yy+=SLICE_LY*dy*SLICE_UP[1];
            zz+=SLICE_LY*dy*SLICE_UP[2];
            
            #if ((GEOMETRY==SPHERICAL) && (SLICE_TYPE==CARTESIAN))
				*rad=pow(xx*xx+yy*yy+zz*zz,0.5);
				*lon=atan2(yy,xx);
				*lat=atan2(zz,pow(xx*xx+yy*yy,0.5));
            #else
				*rad=zz; 
				*lon=xx; 
				*lat=yy; 
			#endif


	#endif
    

}
