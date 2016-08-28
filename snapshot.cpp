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
 

#include <vector>
#include <math.h>   

#include <algorithm>    // std::sort
#include <cstdio>
#include <iostream>
std::string raw_file="raw.dat";

// Must match options.hpp:
//#define QSL_DIM 3
//#define  MAX_RES_BITS 20
#include "options.hpp"

#if QSL_DIM==3
    size_t nx_out=128;
    size_t ny_out=128;
    size_t nz_out=128;
    size_t nout=nx_out*ny_out*nz_out;
#endif
#if QSL_DIM==2
    size_t nx_out=512*8;
    size_t ny_out=512*8;
    size_t nout=nx_out*ny_out;
#endif




//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************





#define  MAX_RES (1<<MAX_RES_BITS) 	 //It should equal 2^MAX_RES_BITS



//#include "options.hpp"
#include "peano.cpp"
struct qsl_struct {
  uint64_t x;
  float qsl;
  size_t i;
  size_t j;
  #if QSL_DIM==3
	size_t k;
  #endif
} ;

struct out_struct {
  float qsl;
} ;

bool comparison (qsl_struct & i,qsl_struct& j) { return (i.x<j.x); }

typedef std::vector< qsl_struct > qsl_type;
typedef std::vector< out_struct > out_type;


int main( void )
{
    using namespace std;
    
    size_t init_size=512*512; // exact value is not important
    
    size_t qsl_size=init_size;
    size_t nin=0;
    


    qsl_type qsl(init_size);
    qsl_type qsl_out(nout);
    out_type arr(nout);
//    vector <double> arr(nout,0.0);
    FILE *file;

    
    //cerr<<"ok0a\n";
    {
        float q,m,bx,by,bz,h,jx,jy,jz,a,v,t;
        uint64_t d;
        file = fopen((raw_file).c_str(), "r");


        size_t i=0;
        float tmp;
        while (!feof(file)) {
				int check=fscanf(file, "%lu%g%g%g%g", &d,&tmp,&tmp,&tmp,&q );
				if ((check!=5)&&(check>0)){
					std::cerr<<"***Problem with raw file at line "<<nin<<'\n';
					exit(0);
			    }
            if (qsl_size==i){
                qsl_size+=init_size;
                qsl.resize(qsl_size);
            }
         #if CALCULATE==QSL
            if (q>1e-6) q=log10(q);       
            else q=-1000;                     
            
            if ((q>100) || (q<-6) || (q!=q)){// Check for junk values.
                q=-1000;
            }
         #endif
            qsl[i].qsl=q;
            qsl[i].x=d;
            nin++;
            i++;
            
        }
        fclose(file);
        qsl.resize(nin);
        std::sort (qsl.begin(), qsl.begin()+nin, comparison);
       }
    
    
    for (uint64_t i=0;i<nx_out;i++)
        for (uint64_t j=0;j<ny_out;j++)
            #if QSL_DIM==3 
                for (uint64_t k=0;k<nz_out;k++)
            #endif
            {
                uint64_t x= (uint64_t) (double(MAX_RES)/double(nx_out)*double(i)+0.5);
                uint64_t y= (uint64_t) (double(MAX_RES)/double(ny_out)*double(j)+0.5);
                #if QSL_DIM==3
                    uint64_t z= (uint64_t) (double(MAX_RES)/double(nz_out)*double(k)+0.5);
                    size_t kk = k+nz_out*(j+ny_out*i);
                    qsl_out[kk].x = peanokey (x,y,z);   
                    qsl_out[kk].k=k;
                #endif
                #if QSL_DIM==2
                    size_t kk = (j+ny_out*i);
                    qsl_out[kk].x = peanokey (x,y);   
                #endif
                qsl_out[kk].qsl=-1000;
                
                qsl_out[kk].i=i;
                qsl_out[kk].j=j;

        }
    std::sort (qsl_out.begin(), qsl_out.end(), comparison);
{
    size_t i=0,k=0;
    size_t k0,k1,i0,i1,t;
    
    uint64_t Hin0,Hin1,Hout0,Hout1;
    
    double d0,d1;
    while (k<nin && i<nout-1){
        Hout0=qsl_out[i].x;
        Hout1=qsl_out[i+1].x;
        
        while (k<nin && qsl[k].x<Hout0) k++;
        //if (k==-1) k=0; 
        k0=k;
        while (k<nin && qsl[k].x<Hout1) k++;
        //std::cerr << k<< " 0 \n";
        if ((k==-1) || (k==0)) k=1; 
        //std::cerr << k<< " 1 \n";
        k--;
        k1=k;
        
        if (k1>=k0){
            for (size_t kk=k0;kk<=k1;kk++){
                if (qsl_out[i].qsl < qsl[kk].qsl){//take the maximum qsl 
                    qsl_out[i].qsl = qsl[kk].qsl;
                   
                }
                 
            }
        }
        else{
            t=k0;k0=k1;k1=t;//swap k0 and k1
            qsl_out[i].qsl=(double(Hout0)-double(qsl[k0].x))/(double(qsl[k1].x)-double(qsl[k0].x)) * (qsl[k1].qsl-qsl[k0].qsl) + qsl[k0].qsl;
                      
        }
        i++;
    }
}
    size_t k;
    for (size_t i=0;i<nout;i++){
        #if QSL_DIM==3
            k=qsl_out[i].k+nz_out*(qsl_out[i].j+ny_out*qsl_out[i].i);
        #endif
        #if QSL_DIM==2
            k=(qsl_out[i].j+ny_out*qsl_out[i].i);
        #endif
        arr[k].qsl=qsl_out[i].qsl;
        
    }
    for (size_t i=0;i<nout;i++){
        cout << arr[i].qsl;   
        cout<< "\n";
    }
    
}
            
