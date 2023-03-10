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

//#define FUSION_MAX_VECTOR_SIZE 12


#include <iostream>
#include <vector>
#include <math.h>   
#include <vexcl/vexcl.hpp>
#include <vexcl/devlist.hpp>
#include <boost/numeric/odeint.hpp>


#include "options.hpp"

int k_bottom;

#include "peano.cpp"

namespace odeint = boost::numeric::odeint;
namespace bnu = boost::numeric::ublas;

//[ vexcl_state_types
typedef vex::vector< double >    vector_type;
typedef vex::vector< cl_double4 >    vector2_type;
typedef vex::vector< cl_ushort4 >    ind2_type;
typedef vex::vector< ushort >    ind_type;
typedef vex::multivector< double, NUM_ODE > state_type;
typedef vex::multivector< double, 3 > b_type;
typedef vex::multivector< ushort, 3 > u_type;
//]

#include <algorithm>    // std::sort
#include <fstream>


int batch_num;
struct qsl_struct {
  uint64_t x;
  bool finishedQ;
  #if (CALCULATE==QSL)
	float qsl;
  #endif 	
  float length;
} ;

bool comparison (const qsl_struct & i,const qsl_struct& j) { return (i.x<j.x); }

typedef std::vector< qsl_struct > qsl_type;

void add_samples_along_hilbert(qsl_type &qsl,size_t *qsl_num_pt);

void initialize_grid( qsl_type &qsl);
void shift_grid( qsl_type &qsl,int64_t shift,size_t qsl_num);
void q_calculation( std::vector< std::vector< std::vector< double >> > &R, vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,qsl_type&qsl,ind2_type &ff_lookup);
void initialize( std::vector< std::vector< std::vector< double >> > &R, vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,qsl_type&qsl,ind2_type &ff_lookup);
#if CALCULATE!=QSL
void reading( vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,ind2_type &ff_lookup,vector2_type &Locals,ind2_type &ff_lookupMid,vector2_type &ffMid);
void integrate_streamlines(std::vector< std::vector< std::vector< double >> > &R,vex::Context &ctx,
			vector2_type &B,vector2_type &ff,ind2_type &ff_lookup,vector2_type &Locals,ind2_type &ff_lookupMid
			,vector2_type &ffMid);
#else
void reading( vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,ind2_type &ff_lookup);
void integrate_streamlines(std::vector< std::vector< std::vector< double >> > &R,vex::Context &ctx,
			vector2_type &B,vector2_type &ff,ind2_type &ff_lookup);
#endif			
void length_calculation( std::vector< std::vector< std::vector< double >> > &R, vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,qsl_type&qsl,ind2_type &ff_lookup);
            



#define DEFS()\
    SIZES() \
    const double z_minimum=solar_radius_sub+ZMIN;

#if CALCULATE>QSL  
std::vector < double > IntegrandA,IntegrandB,IntegrandC,IntegrandD;
std::vector < double > HxMid,HyMid,HzMid;
#endif
std::vector < double > Bx,By,Bz,Hx,Hy,Hz;


#include "interpolation_trilinear.cpp"
#include "find_index.hpp"

DEFS();

#if CALCULATE>QSL             
    #include "find_current.cpp"
    #include "find_hft.cpp"
#endif




const size_t n_box=nx*ny*nz;

#include "lon_lat_rad.cpp"




////
////

////
VEX_FUNCTION(cl_double4,combine,(double, vx)(double, vy)(double,vz ),
    double4 aa={vx,vy,vz,0};
    return aa;
    );
VEX_FUNCTION(cl_double4,combine4,(double, vx)(double, vy)(double,vz )(double,a ),
    double4 aa={vx,vy,vz,a};
    return aa;
    );
VEX_FUNCTION(cl_ushort4,combineu,(ushort, vx)(ushort, vy)(ushort,vz ),
    ushort4 aa={vx,vy,vz,0};
    return aa;
    );
////
VEX_FUNCTION(double,extract0,(cl_double4, v),
    return v.s0;
    );
    
VEX_FUNCTION(double,extract1,(cl_double4, v),
    return v.s1;
    );
    
VEX_FUNCTION(double,extract2,(cl_double4, v),
    return v.s2;
    );

VEX_FUNCTION(double,extract3,(cl_double4, v),
    return v.s3;
    );
    
VEX_FUNCTION(int,extract2u,(cl_ushort4, v),
    int a=v.s2;
    return a;
    );

VEX_FUNCTION(int,extract3u,(cl_ushort4, v),
    int a=v.s3;
    return a;
    );
////
////

#define u1 X(0+3*1)
#define u2 X(1+3*1)
#define u3 X(2+3*1)

#define v1 X(0+3*2)
#define v2 X(1+3*2)
#define v3 X(2+3*2)

double step;

int main( int argc , char **argv )
{
    using namespace std;
    using namespace odeint;
    
    qsl_type qsl(init_size);
    for (size_t i=0;i<init_size;i++)qsl[i].finishedQ=true;

    // setup the opencl context
    vex::Context ctx( vex::Filter::Type(OpenCL_DEVICE_TYPE) && vex::Filter::Position(NGPU) );
    std::cerr << ctx << std::endl;
    
    vector2_type B( ctx.queue() , n_box );
    vector2_type Borig( ctx.queue() , n_box );
    vector2_type ff( ctx.queue() , nmax );
    ind2_type ff_lookup( ctx.queue() , nmax );
    
#if CALCULATE!=QSL
    ind2_type ff_lookupMid( ctx.queue() , nmax );
    vector2_type ffMid( ctx.queue() , nmax );
    vector2_type Locals( ctx.queue() , n_box );
    reading(ctx,B,Borig,ff,ff_lookup,Locals,ff_lookupMid,ffMid);
#else    
    reading(ctx,B,Borig,ff,ff_lookup);
#endif
    cerr<<"Reading successful." << "\n";
    
    double step1;
	step  = (xmax_file-xmin_file)*solar_radius/(((double)nx)*INTEGRATION_STEPS_PER_CELL);
	step1 = (ymax_file-ymin_file)*solar_radius/(((double)ny)*INTEGRATION_STEPS_PER_CELL);
	if (step1<step)step=step1;
	step1 = (zmax_file-zmin_file)/(((double)nz)*INTEGRATION_STEPS_PER_CELL);
	if (step1<step)step=step1;
    cerr<<"Integration step set at "<< step << " Mm\n";
    
    size_t BATCHSIZE=init_size;
    size_t qsl_num=init_size;

    vector< vector< vector< double >> > R;
    R.resize( 2,vector< vector< double >>(NUM_ODE , vector<double>( BATCHSIZE , 0.0 ) ));
    
    initialize_grid(qsl);
    
    
    

    size_t qsl_num_old;
    
    cerr<<"Initialization successful." << "\n";
    
    batch_num=-1;
    double xx,yy,zz;
    
    #ifdef GLOBAL_MODEL
		xmin_file=0.0;
		xmax_file=2.*M_PI;
		ymin_file=-M_PI_2;
		ymax_file=M_PI_2;
        #if CALCULATE!=QSL
		xmin_fileMid=0.0;
		xmax_fileMid=2.*M_PI;
		ymin_fileMid=-M_PI_2;
		ymax_fileMid=M_PI_2;
        #endif
    #endif
    
    while ((BATCHSIZE>0) &&(batch_num<MAX_REFINEMENTS) && (BATCHSIZE<MAX_BATCHSIZE)){
        batch_num++;
        cerr << "Number of field lines to be integrated in this mesh refinement step: " << BATCHSIZE << "\n";
        for (size_t i=0;i<2;i++)
            for (size_t j=0;j<NUM_ODE;j++)
                R[i][j].resize(BATCHSIZE);
        
        {
            size_t k=0,i=0;
            while ((k<BATCHSIZE)&&(i<qsl_num)){
                //cerr<<i<<"  "<<k << "\n";
                while (qsl[i].finishedQ) i++;

                
                hilbert_to_coo(qsl,&xx,&yy,&zz,i);
                #if CALCULATE==QSL
                 if ((xx<=xmax_file) && (xx>=xmin_file) &&
                     (yy<=ymax_file) && (yy>=ymin_file) &&
                     (zz<=zmax_file) && (zz>=zmin_file)) {
				#else
				 if ((xx<=xmax_fileMid) && (xx>=xmin_fileMid) &&
                    (yy<=ymax_fileMid) && (yy>=ymin_fileMid) &&
                    (zz<=zmax_fileMid) && (zz>=zmin_fileMid)) {
				#endif
					R[0][0][k]=xx; 
					R[0][1][k]=yy; 
					R[0][2][k]=zz; 
					R[1][0][k]=xx;
					R[1][1][k]=yy;
					R[1][2][k]=zz;
					k++;
				}
				else{
					qsl[i].finishedQ=true;
					qsl[i].length=0;
					#if CALCULATE==QSL
						qsl[i].qsl=0;
					#endif
					
					#if CALCULATE!=QSL
						cout << qsl[i].x  << "\t" << 0 << "\t" << -1 << "\t" << 0<< "\t" << 0<< "\t" << 0<< "\t" << 0 ;
					#else
						cout << qsl[i].x  << "\t" << -1000 << "\t" << 0 << "\t" << -1 ;
					#endif
                    cout  << "\n";
				}
                i++;
            }
            BATCHSIZE=k;
        }
        
        for (size_t i=0;i<2;i++)
            for (size_t j=0;j<NUM_ODE;j++)
                R[i][j].resize(BATCHSIZE);
        
        
        initialize( R,ctx, B,Borig,ff,qsl,ff_lookup);
        
        
       
			#if (CALCULATE==QSL)
			    integrate_streamlines(R,ctx,B,ff,ff_lookup);
                q_calculation(R,ctx,B,Borig,ff,qsl,ff_lookup);
                cerr<<"Q values calculated successfully." << "\n";
            #else
                integrate_streamlines(R,ctx,B,ff,ff_lookup,Locals,ff_lookupMid,ffMid);
                length_calculation(R,ctx,B,Borig,ff,qsl,ff_lookup);
			    cerr<<"Quantities calculated successfully." << "\n";
            #endif
            
            if (batch_num==MAX_REFINEMENTS) return 0;
        
			BATCHSIZE = 0;
			
			qsl_num_old=qsl_num;
			
			add_samples_along_hilbert(qsl,&qsl_num);
			BATCHSIZE += qsl_num - qsl_num_old;
			
			// Shift grid relative to Hilbert cube, then add samples along hilbert, then shift back.
			// Doing that to fill in more points, which may be missed by the geometry of the hilbert curve.
			// Can skip this if you want to.
			int64_t shift = 3911;//can be anything really. Chosen as an odd number around the more or less random (MAX_RES/256);
			shift_grid(qsl,shift, qsl_num);
			
			qsl_num_old=qsl_num;
			add_samples_along_hilbert(qsl,&qsl_num);
			BATCHSIZE += qsl_num - qsl_num_old;
			
			if (BATCHSIZE>MAX_BATCHSIZE) return 0;
			
			shift = -3911;//shift back
			shift_grid(qsl,shift, qsl_num);
    }
    
    return 0;
}
            
            
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#if CALCULATE==QSL
void q_calculation( std::vector< std::vector< std::vector< double >> > &R, vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,qsl_type&qsl,ind2_type &ff_lookup){
            
            using namespace std;
            size_t batchsize = R[0][0].size();
            size_t n=CHUNKSIZE;
            vector< double > uuF(n);
            vector< double > uuB(n);
            vector< double > vvB(n);
            vector< double > vvF(n);
            vector< double > uvB(n);
            vector< double > uvF(n);
            vector< double > bF(n);
            vector< double > bB(n);
            vector< int > open_fl(n);
            vector< ushort > ind_local(n);
            vector< ushort > ind_k(n);
            
            vector< double > length1(n),length2(n);

            
            b_type b( ctx.queue() , n ); 
            b_type bo( ctx.queue() , n ); 
            vex::vector< double > norm( ctx.queue() , n );
            ind2_type ind(ctx.queue() , n);
            vector2_type d(ctx.queue() , n);
            state_type X(ctx.queue(), n);
            ind_type ind1(ctx.queue(), n);
            
            
            size_t offset=0;
            uint64_t i;
            while (offset<batchsize){
                if (n+offset>batchsize)
                    n=batchsize-offset;
                for (size_t ii = 0; ii < 2; ++ii){
                    for (size_t l = 0; l < NUM_ODE; ++l)
                        vex::copy(R[ii][l].begin() + offset,R[ii][l].begin() + offset +n,X(l).begin());
                    ind = find_index(X(0),X(1),X(2),vex::raw_pointer(ff),vex::raw_pointer(ff_lookup));
                    d = INTERP( X(0),X(1),X(2),vex::raw_pointer(Borig), vex::raw_pointer(ff),ind);
   		            bo(0) = extract0(d);
                    bo(1) = extract1(d);
                    bo(2) = extract2(d);
                    bo(0) = pow(bo(0)*bo(0)+bo(1)*bo(1)+bo(2)*bo(2),0.5)+1.e-30;
                    ind1=extract3u(ind);
                    vex::copy(ind1,ind_local) ;
                    ind1=extract2u(ind);
                    vex::copy(ind1,ind_k) ;
                    d = INTERP( X(0),X(1),X(2),vex::raw_pointer(B), vex::raw_pointer(ff),ind);
                    b(0) = extract0(d);
                    b(1) = extract1(d);
                    b(2) = extract2(d);
                    norm = pow(b(0)*b(0)+b(1)*b(1)+b(2)*b(2),0.5)+1.e-30;
                    b(0) /= norm;
                    b(1) /= norm;
                    b(2) /= norm;
                    #if GEOMETRY==SPHERICAL
                        u1*=X(2)*cos(X(1));
                        u2*=X(2);
                        
                        v1*=X(2)*cos(X(1));
                        v2*=X(2);
                    #endif
                    
                    norm=b(0)*u1+b(1)*u2+b(2)*u3;
                    u1 -= b(0)*norm;
                    u2 -= b(1)*norm;
                    u3 -= b(2)*norm;
                    
                    norm=b(0)*v1+b(1)*v2+b(2)*v3;
                    v1 -= b(0)*norm;
                    v2 -= b(1)*norm;
                    v3 -= b(2)*norm;
                    
            
                    if (ii==0){
						vex::copy(  X(NUM_ODE-1) , length1);
                        b(0) = u1*u1 + u2*u2 + u3*u3;
                        vex::copy(  b(0) , uuF);
                        b(0) = u1*v1 + u2*v2 + u3*v3;
                        vex::copy( b(0) , uvF);
                        b(0) = v1*v1 + v2*v2 + v3*v3;
                        vex::copy(  b(0) , vvF);
                        vex::copy( bo(0) , bF);
                        for (i=0;i<n;i++) {
							bF[i] *= pow(1e9,1.0-((double)ind_local[i]));
							open_fl[i]=0;
							if (ind_k[i]>k_bottom) open_fl[i]=1;
						}
                    }
                    else{
						vex::copy(  X(NUM_ODE-1) , length2);
                        b(0) = u1*u1 + u2*u2 + u3*u3;
                        vex::copy( b(0) , uuB);
                        b(0) = u1*v1 + u2*v2 + u3*v3;
                        vex::copy( b(0) , uvB);
                        b(0) = v1*v1 + v2*v2 + v3*v3;
                        vex::copy( b(0) , vvB);
                        vex::copy( bo(0) , bB);
                        for (i=0;i<n;i++) {
							bB[i] *= pow(1e9, 1.0-((double)ind_local[i]));
							if (ind_k[i]>k_bottom) open_fl[i]=1;
						}
                    }
                
                }
                
                i=0;
                uint64_t k=0;
                double xx,yy,zz;
                
                while (i<n){
                    while (qsl[k].finishedQ) k++;
                    qsl[k].finishedQ=true;
                    qsl[k].length=length1[i]+length2[i];
                    qsl[k].qsl *= ( uuF[i]*vvB[i] + uuB[i]*vvF[i] - 2.0*uvB[i]*uvF[i] ) * (bB[i]*bF[i]) ;
                    hilbert_to_coo(qsl,&xx,&yy,&zz,k);
					cout << qsl[k].x  << "\t" << qsl[k].qsl << "\t" << qsl[k].length << "\t" << open_fl[i] ;
                    cout  << "\n";
		    
                    i++;
                }

                cout.flush();
                offset+=n;
            }
}
#endif
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
void length_calculation( std::vector< std::vector< std::vector< double >> > &R, vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,qsl_type&qsl,ind2_type &ff_lookup){
            
            using namespace std;
            size_t batchsize = R[0][0].size();
            size_t n=CHUNKSIZE;
            vector< ushort > ind_k(n);
            vector< ushort > ind_openF(n),ind_openB(n);
            
            vector< double > length1(n),length2(n);
            vector< double > lengthA1(n),lengthA2(n);
            vector< double > lengthB1(n),lengthB2(n);
            vector< double > lengthC1(n),lengthC2(n);
            vector< double > lengthD1(n),lengthD2(n);
            
            b_type b( ctx.queue() , n ); 
            b_type bo( ctx.queue() , n ); 
            vex::vector< double > norm( ctx.queue() , n );
            ind2_type ind(ctx.queue() , n);
            vector2_type d(ctx.queue() , n);
            state_type X(ctx.queue(), n);
            ind_type ind1(ctx.queue(), n);
            //vector< ushort > ind_local(n);  
            
            size_t offset=0;
            uint64_t i;
            while (offset<batchsize){
                if (n+offset>batchsize)
                    n=batchsize-offset;
                for (size_t ii = 0; ii < 2; ++ii){
                    for (size_t l = 0; l < NUM_ODE; ++l)
                        vex::copy(R[ii][l].begin() + offset,R[ii][l].begin() + offset +n,X(l).begin());
                    ind = find_index(X(0),X(1),X(2),vex::raw_pointer(ff),vex::raw_pointer(ff_lookup));
                    ind1=extract2u(ind);
                    vex::copy(ind1,ind_k) ;

                    if (ii==0){
						vex::copy(  X(7) , lengthD1);
						vex::copy(  X(6) , lengthC1);
						vex::copy(  X(5) , lengthB1);
						vex::copy(  X(4) , lengthA1);
						vex::copy(  X(3) , length1);
                        for (i=0;i<n;i++) {
							ind_openF[i]=0;
								if (ind_k[i]>k_bottom) ind_openF[i]=1;
						}
                    }
                    else{
						vex::copy(  X(7) , lengthD2);
						vex::copy(  X(6) , lengthC2);
						vex::copy(  X(5) , lengthB2);
						vex::copy(  X(4) , lengthA2);
						vex::copy(  X(3) , length2);
                        for (i=0;i<n;i++) {
							ind_openB[i]=0;
							if (ind_k[i]>k_bottom) ind_openB[i]=1;
						}
                    }
                
                }
                
                i=0;
                uint64_t k=0;
                double xx,yy,zz;
                
                while (i<n){
                    while (qsl[k].finishedQ) k++;
                    qsl[k].finishedQ=true;
                      double a,b,c,d;
					  d=(lengthD1[i]+lengthD2[i]);
					  c=(lengthC1[i]+lengthC2[i]);
					  b=(lengthB1[i]+lengthB2[i]);
					  a=(lengthA1[i]+lengthA2[i]);
					  qsl[k].length=(length1[i]+length2[i]);
					  
					  
					  //*((double)ind_openB[i]*ind_openF[i]);
					  hilbert_to_coo(qsl,&xx,&yy,&zz,k);
					  //<< "\t"<< xx<<"\t"<< yy<<"\t"<< zz
					  cout << qsl[k].x  << "\t" << qsl[k].length << "\t" << 1-(1-ind_openB[i])*(1-ind_openF[i]) << "\t" << a<< "\t" << b<< "\t" << c<< "\t" << d ; // hilbert, x,y,z,length          
                    cout  << "\n";
		    
                    i++;
                }

                cout.flush();
                offset+=n;
            }
}

        


void initialize( std::vector< std::vector< std::vector< double >> > &R, vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,qsl_type&qsl,ind2_type &ff_lookup){
        using namespace std;
        size_t n=CHUNKSIZE;
        state_type X(ctx.queue(), n);
        b_type b( ctx.queue() , n ); 
        vex::vector< double > norm( ctx.queue() , n );
        std::vector< double > BB(n);
        
        ind2_type ind(ctx.queue() , n);
        vector2_type d(ctx.queue() , n);
        size_t batchsize = R[0][0].size();
        ind_type ind1(ctx.queue(), n);
        vector< ushort > ind_local(n);  
        
        size_t offset=0;
        while (offset<batchsize){
            if (n+offset>batchsize)
                n=batchsize-offset;
            vex::copy(R[0][0].begin() + offset,R[0][0].begin()+n + offset,X(0).begin());
            vex::copy(R[0][1].begin() + offset,R[0][1].begin()+n + offset,X(1).begin());
            vex::copy(R[0][2].begin() + offset,R[0][2].begin()+n + offset,X(2).begin());

#if CALCULATE==QSL
            ind = find_index(X(0),X(1),X(2),vex::raw_pointer(ff),vex::raw_pointer(ff_lookup));
            ind1=extract3u(ind);
            vex::copy(ind1,ind_local) ;
	        d =  INTERP( X(0),X(1),X(2), vex::raw_pointer(Borig), vex::raw_pointer(ff),ind);
            b(0) = extract0(d) ; 
            b(1) = extract1(d) ; 
            b(2) = extract2(d) ; 
            
            b(0)=b(0)*b(0)+b(1)*b(1)+b(2)*b(2);
            vex::copy(b(0).begin(),b(0).begin()+n,BB.begin()); ////
            {
                uint64_t k=0,i=0;
                for (i=0;i<offset;i++){
                    while (qsl[k].finishedQ) k++;
                    k++;
                }
                i=0;
                while (i<n){
                    while (qsl[k].finishedQ) k++;
                    qsl[k].qsl=1./BB[i]/pow(1e18,1.0-((double)ind_local[i]));
                    qsl[k].length=0.0;
                    i++;
                    k++;
                }
            }
	        d =  INTERP( X(0),X(1),X(2), vex::raw_pointer(B), vex::raw_pointer(ff),ind);
            //setup IC for perturbations using Gram-Schmidt
            
            b(0) = extract0(d) ; 
            b(1) = extract1(d) ; 
            b(2) = extract2(d) ; 
            

            norm = pow(b(0)*b(0)+b(1)*b(1)+b(2)*b(2),0.5)+1.e-30;
            b(0) /= norm;
            b(1) /= norm;
            b(2) /= norm;
            
            // now x,y,z has all the b vector
            X(0+3*2) = b(1);
            X(1+3*2) = -b(0);
            X(2+3*2) = 0.0;
            norm =  pow(X(0+3*2)* X(0+3*2)+ X(1+3*2)*X(1+3*2),0.5)+1.e-30;
            X(0+3*2) /= norm; // A = y   normalized
            X(1+3*2) /= norm; // B = -x  normalized
            // now 3*2 has all the tildas
            X(0+3*1) = -X(1+3*2)*b(2); // x=-B*z
            X(1+3*1) =  X(0+3*2)*b(2);  // y=A*z
            X(2+3*1) = - X(0+3*2)*b(1) + X(1+3*2)*b(0); // C = -A*y + B*x
            
            #if GEOMETRY==SPHERICAL
				X(0+3*1)/=X(2)*cos(X(1));
				X(1+3*1)/=X(2);
				X(0+3*2)/=X(2)*cos(X(1));
				X(1+3*2)/=X(2);
            #endif
            X(9)=0;
#else
            {
                uint64_t k=0,i=0;
                for (i=0;i<offset;i++){
                    while (qsl[k].finishedQ) k++;
                    k++;
                }
                i=0;
                while (i<n){
                    while (qsl[k].finishedQ) k++;
                    qsl[k].length=0.0;
                    i++;
                    k++;
                }
            }
			for (size_t l = 3; l < NUM_ODE; ++l)
				X(l)=0;
#endif
            for (size_t l = 3; l < NUM_ODE; ++l){
                vex::copy(X(l).begin(),X(l).begin()+n, R[0][l].begin() + offset);
                vex::copy(X(l).begin(),X(l).begin()+n, R[1][l].begin() + offset);
            }
            offset+=n;
            
        }
    }



#include <cstdio>
#if CALCULATE!=QSL
void reading( vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,ind2_type &ff_lookup,vector2_type &Locals,ind2_type &ff_lookupMid,vector2_type &ffMid){
#else
void reading( vex::Context &ctx, 
                    vector2_type &B,vector2_type &Borig,vector2_type &ff,ind2_type &ff_lookup){
#endif
    using namespace std;
    b_type bbf( ctx.queue() , nmax ); 
    u_type bbfu( ctx.queue() , nmax ); 
    //vector2_type ff(ctx.queue(), nmax);
	if ((nx>=nmax) || (ny>=nmax) || (nz>=nmax)) {
		std::cerr<<"***Increase nmax in post-defs.hpp. Currently, nmax="<<nmax <<", while nx,ny,nz="<<nx<<' '<<ny<<' '<<nz<<'\n';
		exit(0);
	}
    
    FILE *filefx;
    FILE *filefy;
    FILE *filefz;
    FILE *filex;
    FILE *filey;
    FILE *filez;
    
      

    
    Hx.resize(  nmax, 0.0 );
    Hy.resize(  nmax, 0.0 );
    Hz.resize(  nmax, 0.0 );
    
    vector< ushort > fx_lookup(nmax);
    vector< ushort > fy_lookup(nmax);
    vector< ushort > fz_lookup(nmax);
#if CALCULATE!=QSL
    HxMid.resize(  nmax, 0.0 );
    HyMid.resize(  nmax, 0.0 );
    HzMid.resize(  nmax, 0.0 );
    
    vector< ushort > fx_lookupMid(nmax);
    vector< ushort > fy_lookupMid(nmax);
    vector< ushort > fz_lookupMid(nmax);
#endif

    double ox,oy,oz,o;
    int check,i; 
    filefx = fopen((in_dir+"xs0"+in_filename+".dat").c_str(), "r");
    if (filefx==0){
		std::cerr<<"***Missing xs0 file"<<'\n';
		exit(0);
	}
    check=1;i=0;
    while (check==1) {
        check=fscanf(filefx, "%lf", &o );
        if (check==1)Hx[i]=o*to_radians;
        i++;
    }
    if (nx!=i-1) {
		std::cerr<<"***MISMATCHED SIZES for xs0"<<'\n';
		exit(0);
	}
    xmax_file=Hx[nx-1];
	xmin_file=Hx[0];
#if QSL_DIM==3
	#ifndef GLOBAL_MODEL
    if (xmin_file/to_radians>XMIN) {
		std::cerr<<    "***XMIN OFF: "<<xmin_file/to_radians<<" "    <<XMIN<<'\n';
		exit(0);
	}
    if (xmax_file/to_radians<XMAX) {
		std::cerr<<"***XMAX OFF: "<<xmax_file/to_radians<<" "<<XMAX<<'\n';
		exit(0);
	}
	#else
	if (xmin_file/to_radians<0) {
		std::cerr<<    "***XMIN OFF: "<<xmin_file/to_radians<<" "    <<0<<'\n';
		exit(0);
	}
    if (xmax_file/to_radians>360) {
		std::cerr<<"***XMAX OFF: "<<xmax_file/to_radians<<" "<<360<<'\n';
		exit(0);
	}
	#endif
	
	
	
	
	
#endif
    fclose(filefx);
    vex::copy( Hx,bbf(0));
    
    
    
    
    
    filefy = fopen((in_dir+"ys0"+in_filename+".dat").c_str(), "r");
    if (filefy==0){
		std::cerr<<"***Missing ys0 file"<<'\n';
		exit(0);
	}
    check=1;i=0;
    while (check==1) {
        check=fscanf(filefy, "%lf", &o );
        if (check==1)Hy[i]=o*to_radians;
        i++;
    }
    if (ny!=i-1) {
		std::cerr<<"***MISMATCHED SIZES for ys0"<<'\n';
		exit(0);
	}
    ymax_file=Hy[ny-1];
	ymin_file=Hy[0];
#if QSL_DIM==3
	#ifndef GLOBAL_MODEL
    if (ymin_file/to_radians>YMIN) {
		std::cerr<<"***YMIN OFF: "<<ymin_file/to_radians<<" "<<YMIN<<'\n';
		exit(0);
	}
    if (ymax_file/to_radians<YMAX) {
		std::cerr<<"***YMAX OFF: "<<ymax_file/to_radians<<" "<<YMAX<<'\n';
		exit(0);
	}
	#else
	if (ymin_file/to_radians<-90) {
		std::cerr<<"***YMIN OFF: "<<ymin_file/to_radians<<" "<<-90<<'\n';
		exit(0);
	}
    if (ymax_file/to_radians>90) {
		std::cerr<<"***YMAX OFF: "<<ymax_file/to_radians<<" "<<90<<'\n';
		exit(0);
	}
	#endif
#endif
    fclose(filefy);
    vex::copy( Hy,bbf(1));
    
    filefz = fopen((in_dir+"zs0"+in_filename+".dat").c_str(), "r");
    if (filefy==0){
		std::cerr<<"***Missing zs0 file"<<'\n';
		exit(0);
	}
    check=1;i=0;k_bottom=0;
    while (check==1) {
        check=fscanf(filefz, "%lf", &o );
        if (check==1){
			Hz[i]=o*solar_radius;
			if (Hz[i]-solar_radius_sub<=ZMIN) k_bottom=i;
		}
        i++;
    }
    k_bottom++;
    if (nz!=i-1) {
		std::cerr<<"***MISMATCHED SIZES for zs0"<<'\n';
		exit(0);
	}
	zmax_file=Hz[nz-1];
	zmin_file=Hz[0];
#if QSL_DIM==3
    if (zmin_file-solar_radius_sub>ZMIN) {
		std::cerr<<    "***ZMIN OFF: "<<zmin_file-solar_radius_sub<<" "    <<ZMIN<<'\n';
		exit(0);
	}
    if (zmax_file-solar_radius_sub<ZMAX) {
		std::cerr<<"***ZMAX OFF: "<<    zmax_file-solar_radius_sub<<" "<<ZMAX<<'\n';
        exit(0);
	}
#endif
    fclose(filefz);
    vex::copy( Hz,bbf(2));
    
    ff = combine(bbf(0),bbf(1),bbf(2));
#if CALCULATE!=QSL
    for (size_t k = 0; k < nz-1; ++k)  
		HzMid[k]=(Hz[k]+Hz[k+1])/2.;
	for (size_t k = 0; k < ny-1; ++k)  
		HyMid[k]=(Hy[k]+Hy[k+1])/2.;
    for (size_t k = 0; k < nx-1; ++k)  
		HxMid[k]=(Hx[k]+Hx[k+1])/2.;
	zmax_fileMid=HzMid[nz-2];
	zmin_fileMid=HzMid[0];
	ymax_fileMid=HyMid[ny-2];
	ymin_fileMid=HyMid[0];
	xmax_fileMid=HxMid[nx-2];
	xmin_fileMid=HxMid[0];
	
	vex::copy( HxMid,bbf(0));
	vex::copy( HyMid,bbf(1));
	vex::copy( HzMid,bbf(2));
    
    ffMid = combine(bbf(0),bbf(1),bbf(2));
#endif
///// CREATE LOOKUP TABLES
    ushort kx=0,ky=0,kz=0;
    ushort kxM=0,kyM=0,kzM=0;
    for (size_t i = 0; i < nmax; i++){
		double ix,iy,iz;
		ix=(xmax_file-xmin_file)*((double)i)/((double)nmax)+xmin_file;
		iy=(ymax_file-ymin_file)*((double)i)/((double)nmax)+ymin_file;
		iz=(zmax_file-zmin_file)*((double)i)/((double)nmax)+zmin_file;
		while ((kx<nx) && (ix>Hx[kx])) kx++;
		fx_lookup[i]=kx;
		while ((ky<ny) && (iy>Hy[ky])) ky++;
		fy_lookup[i]=ky;
		while ((kz<nz) && (iz>Hz[kz])) kz++;
		fz_lookup[i]=kz;
		
#if CALCULATE!=QSL
		ix=(xmax_fileMid-xmin_fileMid)*((double)i)/((double)nmax)+xmin_fileMid;
		iy=(ymax_fileMid-ymin_fileMid)*((double)i)/((double)nmax)+ymin_fileMid;
		iz=(zmax_fileMid-zmin_fileMid)*((double)i)/((double)nmax)+zmin_fileMid;
		while ((kxM<nx-1) && (ix>HxMid[kxM])) kxM++;
		fx_lookupMid[i]=kxM;
		while ((kyM<ny-1) && (iy>HyMid[kyM])) kyM++;
		fy_lookupMid[i]=kyM;
		while ((kzM<nz-1) && (iz>HzMid[kzM])) kzM++;
		fz_lookupMid[i]=kzM;
#endif
	}
    
    vex::copy( fx_lookup,bbfu(0));
    vex::copy( fy_lookup,bbfu(1));
    vex::copy( fz_lookup,bbfu(2));
    
    ff_lookup = combineu(bbfu(0),bbfu(1),bbfu(2));
#if CALCULATE!=QSL
    vex::copy( fx_lookupMid,bbfu(0));
    vex::copy( fy_lookupMid,bbfu(1));
    vex::copy( fz_lookupMid,bbfu(2));
    
    ff_lookupMid = combineu(bbfu(0),bbfu(1),bbfu(2));
#endif    
    
///////
    // B( ctx.queue() , n_box );   
    b_type bb( ctx.queue() , n_box ); 
///////
    vector< double > datax(n_box);
    vector< double > datay(n_box);
    vector< double > dataz(n_box);

    Bx.resize(  n_box, 0.0 );
    By.resize(  n_box, 0.0 );
    Bz.resize(  n_box, 0.0 );


    filex = fopen((in_dir+"bx0"+in_filename+".dat").c_str(), "r");
    filey = fopen((in_dir+"by0"+in_filename+".dat").c_str(), "r");
    filez = fopen((in_dir+"bz0"+in_filename+".dat").c_str(), "r"); 
    if (filex==0){
		std::cerr<<"*** Missing bx0 file."<<'\n';
        exit(0);
	}
	if (filey==0){
		std::cerr<<"*** Missing by0 file."<<'\n';
        exit(0);
	}
    if (filez==0){
		std::cerr<<"*** Missing bz0 file."<<'\n';
        exit(0);
	}
    check=1;
    for (size_t k = 0; k < nz; ++k)  
        for (size_t j= 0; j < ny; ++j) 
            for (size_t i = 0; i < nx; ++i)  {
				
                check=fscanf(filex, "%lf", &ox )*check;
                check=fscanf(filey, "%lf", &oy )*check;
                check=fscanf(filez, "%lf", &oz )*check;
                if (check==1){
					o=pow(ox*ox+oy*oy+oz*oz,0.5)+1.e-30;
					datax[i+ nx*(j+ny*k)]=ox   /o; //*10.0;
					datay[i+ nx*(j+ny*k)]=oy   /o; //*10.0;
					dataz[i+ nx*(j+ny*k)]=oz   /o; //*10.0;
					Bx[i+nx*(j+ny*k)]=ox;//   /o; //*10.0;
					By[i+nx*(j+ny*k)]=oy;//   /o; //*10.0;
					Bz[i+nx*(j+ny*k)]=oz;//   /o; //*10.0;
				}
            }
    if (check!=1){
		std::cerr<<"*** Problem with B field files."<<'\n';
        exit(0);
	}
    fclose(filex);
    fclose(filey);
    fclose(filez);
    
    
    vex::copy( datax,bb(0) );
    vex::copy( datay,bb(1) );
    vex::copy( dataz,bb(2) );
    

#if CALCULATE>QSL
    IntegrandA.resize(n_box,0.);
    IntegrandB.resize(n_box,0.);
    IntegrandC.resize(n_box,0.);
    IntegrandD.resize(n_box,0.);
    {
    vector_type aaa( ctx.queue() , n_box ); 
    vector_type bbb( ctx.queue() , n_box ); 
    vector_type ccc( ctx.queue() , n_box ); 
    vector_type ddd( ctx.queue() , n_box ); 
	
    find_current();
    find_hft();
    vex::copy( IntegrandA,aaa );
    vex::copy( IntegrandB,bbb );
    vex::copy( IntegrandC,ccc );
    vex::copy( IntegrandD,ddd );
    Locals = combine4(aaa,bbb,ccc,ddd);
	}
#endif
    B = combine(bb(0),bb(1),bb(2));


    
    vex::copy( Bx,bb(0) );
    vex::copy( By,bb(1) );
    vex::copy( Bz,bb(2) );
    Borig = combine(bb(0),bb(1),bb(2));


    
}


#if CALCULATE!=QSL
void integrate_streamlines(std::vector< std::vector< std::vector< double >> > &R,vex::Context &ctx,
			vector2_type &B,vector2_type &ff,ind2_type &ff_lookup,vector2_type &Locals,ind2_type &ff_lookupMid
			,vector2_type &ffMid){
#else
void integrate_streamlines(std::vector< std::vector< std::vector< double >> > &R,vex::Context &ctx,
			vector2_type &B,vector2_type &ff,ind2_type &ff_lookup){
#endif
    using namespace std;
    
    size_t counter,start;
    
    vector< vector< double >> tmp;
    tmp.resize( NUM_ODE , vector<double>( CHUNKSIZE , 0.0 ) );
    size_t BATCHSIZE=  R[0][0].size();
    
    
    vector<size_t> Rind(CHUNKSIZE);
    vector<size_t> R_not_dispatched_yet(BATCHSIZE);
    vector<double> field_line_length(CHUNKSIZE);
    
    size_t n;
    for (size_t ii=0; ii<2;ii++){// loop over forward and backward integration along streamline
		    
            if (ii==1){
				B*=-1.0;//reverse integration direction. i.e. follow streamline backwards
				cerr << "Beginning BACKWARD integration along field lines ...\n";
			}
			else
				cerr << "Beginning FORWARD integration along field lines ...\n";
            fill(R_not_dispatched_yet.begin(),R_not_dispatched_yet.end(),1); // 1= not started
            counter=0;
            n=CHUNKSIZE;
            start=0;
            while (BATCHSIZE > counter){
                cerr<< "# of computed field lines = "<< counter <<" out of "<<BATCHSIZE <<" in mesh refinement: "<<batch_num<<"\n";
                if (BATCHSIZE < (CHUNKSIZE+ counter))
                    n =  BATCHSIZE -counter; // should not try to send more vectors to GPU than the available number in batchsize
                state_type X(ctx.queue(), n); // recreate X every time -- long story ...
                state_type dxdt(ctx.queue(), n);
                vex::vector< int > cQ(ctx.queue(), n);
                vector< int >convergedQ(n);
                
                {
                    size_t k=0;
                    size_t i=start; 
                    // Some of the stream lines have not terminated at box boundaries yet; 
                    // so need to resend them to GPU. 
                    // Those are stored in Rind up to i=start.
                    // Refill Rind above i=start with new IC's that have not been 
                    // dispatched yet.
                    while((k<BATCHSIZE) && (i<n)){
                        if (R_not_dispatched_yet[k]==1){//dispatch the k-th IC
                            Rind[i]=k;
                            field_line_length[i]=0.0; 
                            R_not_dispatched_yet[k]=0;
                            i++;
                        }
                        k++;
                    }
                }
                
                for (size_t l=0; l<NUM_ODE;l++){
                        for (size_t i=start; i<n;i++)
                            tmp[l][i]=R[ii][l][Rind[i]];
                        vex::copy(tmp[l].begin() ,tmp[l].begin()+n,X(l).begin());
                }
		
	            /////////////////////////////////
				/////////////////////////////////
				//////Euler integration//////////
				/////////////////////////////////
				/////////////////////////////////
				/////////////////////////////////
				{	
					
					VEX_CONSTANT(dt, step);
					int steps=(int)(INTEGRATION_RANGE/step);
					if (steps<=0){
						std::cerr<<"***Integration range and timestep result in number of steps<=0. steps = "<<   steps << '\n';
						exit(0);
					}
					for (size_t y=0; y<steps;y++){
						#if (GEOMETRY==SPHERICAL)
							//fix periodicities
							X(0)=fmod(X(0)-xmin_file+(floor(fabs(X(1))*M_2_PI)+2.)*M_PI,2.*M_PI)+xmin_file;
							X(1)=M_PI_2-fabs(fabs(X(1)+M_PI_2)-M_PI);
						#else
							    #ifdef PERIODIC_XY
								X(0)=fmod(X(0)-xmin_file+(xmax_file-xmin_file),xmax_file-xmin_file)+xmin_file;
								X(1)=fmod(X(1)-ymin_file+(ymax_file-ymin_file),ymax_file-ymin_file)+ymin_file;
								
								#endif
						#endif
						
						
						
						ind2_type ind = find_index(X(0),X(1),X(2),vex::raw_pointer(ff),vex::raw_pointer(ff_lookup));
                        vector2_type d =  INTERP( X(0),X(1),X(2), vex::raw_pointer(B), vex::raw_pointer(ff),ind);
		                #if (GEOMETRY==SPHERICAL)
                            dxdt(0) = extract0(d)/(X(2)*cos(X(1)));  
                            dxdt(1) = extract1(d)/X(2); 
                        #else
                            dxdt(0) = extract0(d);
                            dxdt(1) = extract1(d);
                        #endif
						dxdt(2) = extract2(d);  
						  
#if (CALCULATE>QSL)	
						dxdt(3) = extract3(d);
						dxdt(0)*=dxdt(3);
						dxdt(1)*=dxdt(3);
						dxdt(2)*=dxdt(3);
						ind = find_indexMid(X(0),X(1),X(2),vex::raw_pointer(ffMid),vex::raw_pointer(ff_lookupMid));
                        d =  INTERPMID( X(0),X(1),X(2), vex::raw_pointer(Locals), vex::raw_pointer(ffMid),ind);
						dxdt(4) = extract0(d)*dxdt(3);  
						dxdt(5) = extract1(d)*dxdt(3);  
						dxdt(6) = extract2(d)*dxdt(3);  
						dxdt(7) = extract3(d)*dxdt(3); 
#endif		  
#if (CALCULATE==QSL)	
						dxdt(9) = extract3(d);  
						dxdt(0)*=dxdt(9);
						dxdt(1)*=dxdt(9);
						dxdt(2)*=dxdt(9);
						d = INTERP_DIFF(X(0),X(1),X(2),X(0+3*1),X(1+3*1),X(2+3*1), vex::raw_pointer(B), vex::raw_pointer(ff),ind);	
						#if (GEOMETRY==SPHERICAL)
                            dxdt(0+3*1) = extract0(d)*dxdt(9)/(X(2)*cos(X(1)));  
                            dxdt(1+3*1) = extract1(d)*dxdt(9)/X(2); 
                            dxdt(1+3*1) -= X(2+3*1)/X(2)*dxdt(1); 
                            dxdt(0+3*1) -= (X(2+3*1)/X(2)-X(1+3*1)*tan(X(1)))*dxdt(0);
                        #else
                            dxdt(0+3*1) = extract0(d)*dxdt(9);
                            dxdt(1+3*1) = extract1(d)*dxdt(9);
                        #endif
						dxdt(2+3*1) = extract2(d)*dxdt(9); 
                        
                        d = INTERP_DIFF(X(0),X(1),X(2),X(0+3*2),X(1+3*2),X(2+3*2), vex::raw_pointer(B), vex::raw_pointer(ff),ind);
						#if (GEOMETRY==SPHERICAL)
                            dxdt(0+3*2) = extract0(d)*dxdt(9)/(X(2)*cos(X(1)));  
                            dxdt(1+3*2) = extract1(d)*dxdt(9)/X(2); 
                            dxdt(1+3*2) -= X(2+3*2)/X(2)*dxdt(1);
                            dxdt(0+3*2) -= (X(2+3*2)/X(2)-X(1+3*2)*tan(X(1)))*dxdt(0);
                        #else
                            dxdt(0+3*2) = extract0(d)*dxdt(9);
                            dxdt(1+3*2) = extract1(d)*dxdt(9);
                        #endif
						dxdt(2+3*2) = extract2(d)*dxdt(9);
#endif
						for (uint l=0;l<NUM_ODE;l++)
							X(l)+=dxdt(l)*dt();
					}
					
					
					
				}
				//////////////////////////////////////
				///////////////////////////////////////
				//////////////////////////////////////////
				//////////////////////////////////////////
				///////////////////////////////////////////
				/////////////////////////////////

				#if (GEOMETRY==SPHERICAL)
					//fix periodicities
					X(0)=fmod(X(0)-xmin_file+(floor(fabs(X(1))*M_2_PI)+2.)*M_PI,2.*M_PI)+xmin_file;
					X(1)=M_PI_2-fabs(fabs(X(1)+M_PI_2)-M_PI);
					
				#else
					    #ifdef PERIODIC_XY
								X(0)=fmod(X(0)-xmin_file+(xmax_file-xmin_file),xmax_file-xmin_file)+xmin_file;
								X(1)=fmod(X(1)-ymin_file+(ymax_file-ymin_file),ymax_file-ymin_file)+ymin_file;
						#endif
				#endif
		
		
				ind2_type ind(ctx.queue() , n);
                ind=find_index(X(0),X(1),X(2),vex::raw_pointer(ff),vex::raw_pointer(ff_lookup));
                cQ =   extract3u(ind);
                /////////
                //check for small velocity
                vector2_type d(ctx.queue() , n);
                std::vector< double > BB(n);
				b_type b( ctx.queue() , n ); 

				d =  INTERP( X(0),X(1),X(2), vex::raw_pointer(B), vex::raw_pointer(ff),ind);

				b(0) = extract0(d) ; 
				b(1) = extract1(d) ; 
				b(2) = extract2(d) ; 
				b(0)=b(0)*b(0)+b(1)*b(1)+b(2)*b(2);
				vex::copy(b(0).begin(),b(0).begin()+n,BB.begin()); 
                /////////
                
                vex::copy(cQ,convergedQ);
                
                
                for (size_t l = 0; l < NUM_ODE; ++l){
                    start=0;
                    vex::copy(X(l).begin(),X(l).begin()+n,tmp[l].begin());
                    for (size_t i = 0; i < n; ++i) {
						    //cout <<i<<"\t"<< BB[i] << "\n";
                            #ifndef LOCAL_Q 
                                if ((convergedQ[i]==0) ||(BB[i] < 1e-12) || (!(BB[i]==BB[i])) || field_line_length[i] >=MAX_HALF_LENGTH_FIELD_LINE ){ 
									    // zero stands for converged
										// sometimes nan's appear in BB (seems to be a problem only for single precision. I think it's due to not enough precision when subtracting 1e-6 from boundaries ....)
										// get rid of them by checking for ineq
                            #endif
                                R[ii][l][Rind[i]]=tmp[l][i];
                                if(l==0){
                                    counter++;
                                }
                            #ifndef LOCAL_Q 
                                }
                                else{// move all non-terminated field lines to index i=start
                                    tmp[l][start]=tmp[l][i];
                                    if(l==(NUM_ODE-1)){
                                        Rind[start]=Rind[i];
                                        field_line_length[start]=field_line_length[i]+INTEGRATION_RANGE;
                                    }
                                    start++;
                                }
                            #endif
                    }
                }
                

                //cout<<"ok "<<sum(cQ) << " "<< start <<"  "<<counter<<" \n";
                

                

            }
        }
        B*=-1.0;//restore sign of B
}




//////////////////////////////////
void add_samples_along_hilbert(qsl_type &qsl,size_t *qsl_num_pt){
        
        uint64_t x,x1,k1,d1,k,dx; //,m0,m1;
        size_t qsl_num = qsl_num_pt[0];
        size_t qsl_size=qsl.size();
        std::cerr<<"Starting sort..." << "\n";
        std::sort (qsl.begin(), qsl.begin()+qsl_num, comparison);
        std::cerr<<"... done sorting." << "\n";
        bool refine;
        double qC,qP;
        
        k=0;
        while (k<qsl_num){
            while (!qsl[k].finishedQ)k++;
            k1=k+1;
            while (k1!=qsl_num && !qsl[k1].finishedQ) k1=k1+1;

            if (k1!=qsl_num){
				qC=qsl[k].length;
				qP=qsl[k1].length;
				refine =1;
				#if CALCULATE==QSL
					refine = (qsl[k].qsl>2.) && (qsl[k1].qsl>2.);//NO LONGER VALID as q is not overriden by FL openness, but check is useful for sane numbers.
                #endif
                refine = refine && (fabs(qC-qP) > LENGTH_JUMP_REFINEMENT_THRESHOLD);
                refine = refine && (qC>0) && (qP>0);
                if (refine) {
                    x  = qsl[k].x;
                    x1 = qsl[k1].x;
                    if (x1 >x) {
						dx=x1-x;
						d1=x+dx/2;
					}
                    else{
						dx=x-x1;
						d1=x1+dx/2;
					}
                    if ((d1!=x) && (d1!=x1) && (k1==k+1)) {
                        if (qsl_size==qsl_num){
                            qsl_size+=init_size;
                            qsl.resize(qsl_size);
                            for (size_t i=qsl_num;i<qsl_size;i++)qsl[i].finishedQ=true;
                        }
                        qsl[qsl_num].x = d1;
                        qsl[qsl_num].finishedQ=false;
                        qsl_num++;
                    }
                }
            }
            k++;
        }

        qsl_num_pt[0] = qsl_num;
    }

/////////////////

void initialize_grid( qsl_type &qsl){
    
    size_t iq=0;
    
    for (uint64_t i=0;i<nx_init;i++)
        for (uint64_t j=0;j<ny_init;j++)
            #if QSL_DIM==3
                for (uint64_t k=0;k<nz_init;k++)
            #endif
        {
            uint64_t x= (uint64_t) (double(MAX_RES)/double(nx_init-1)*double(i)+0.5);
            uint64_t y= (uint64_t) (double(MAX_RES)/double(ny_init-1)*double(j)+0.5);
            if (x==MAX_RES) x=MAX_RES-1;
            if (y==MAX_RES) y=MAX_RES-1;
            //size_t k = zz+nz_init*(j+ny_init*i);
            #if QSL_DIM==3 
                uint64_t z= (uint64_t) (double(MAX_RES)/double(nz_init-1)*double(k)+0.5);
                if (z==MAX_RES) z=MAX_RES-1;
                qsl[iq].x = peanokey (x,y,z);
            #endif
            #if QSL_DIM==2
                qsl[iq].x = peanokey (x,y);
            #endif
            qsl[iq].finishedQ = false;
            iq++;
        }
}



void shift_grid( qsl_type &qsl,int64_t shift,size_t qsl_num){
    
    uint64_t x,y,z=0; // z is not used in 3d
    for (uint64_t i=0;i<qsl_num;i++){
        #if QSL_DIM==3
			point( qsl[i].x , &x ,&y,&z);
        #endif
        #if QSL_DIM==2
			point( qsl[i].x , &x ,&y);
        #endif
            x+=MAX_RES;
            y+=MAX_RES;
            z+=MAX_RES;
            
            if (shift>0){
                y+=shift;
                x+=shift;
                z+=shift;
            }
            else {
                z-=(uint64_t)(-shift);
                y-=(uint64_t)(-shift);
                x-=(uint64_t)(-shift);
            }
            x%=MAX_RES;
            y%=MAX_RES;
            z%=MAX_RES;
            
        #if QSL_DIM==3
            qsl[i].x = peanokey(x,y,z);
	    #endif
        #if QSL_DIM==2
            qsl[i].x = peanokey(x,y);
	    #endif
        }
}
