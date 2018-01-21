/*
	This file is part of QSL Squasher. 
	Copyright (C) 2014-2018  Svetlin Tassev
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

struct out_struct {
  float qsl;
  int gap;
} ;

typedef std::vector< out_struct > out_type;


int main( void )
{
    using namespace std;
    size_t nin=0;
    out_type arr(nout);
    FILE *file;
	for (size_t kk=0;kk<nout;kk++){
		arr[kk].qsl=-1000;
		arr[kk].gap=-1;
	}
    //cerr<<"ok0a\n";
    {
        float q,m,bx,by,bz,h,a,v,t;
        uint64_t d,k;
        int ofl;
        uint64_t xx,yy,zz;
        uint64_t xx_out,yy_out,zz_out;
        file = fopen((raw_file).c_str(), "r");


        float tmp;
        while (!feof(file)) {
                nin+=1;
				int check=fscanf(file, "%lu%g%g%g%g", &d,&tmp,&tmp,&tmp,&q );
				if ((check!=5)&&(check>0)){
					std::cerr<<"***Problem with raw file at line "<<nin<<'\n';
					exit(0);
			    }
         #if CALCULATE==QSL
            if (q>1e-6) q=log10(q);       
            else q=-1000;                     
            
            if ((q>100) || (q<-6) || (q!=q)){// Check for junk values.
                q=-1000;
            }
         #endif
         
            #if (QSL_DIM == 3)
				point( d, &xx, &yy, &zz);
				zz_out=(uint64_t) (double(zz)/double(MAX_RES-1)*double(nz_out-1)+0.5);
				if (zz_out>nz_out-1)zz_out=nz_out-1;
			#else
				point( d, &xx, &yy);
			#endif
           
            xx_out=(uint64_t) (double(xx)/double(MAX_RES-1)*double(nx_out-1)+0.5);
            yy_out=(uint64_t) (double(yy)/double(MAX_RES-1)*double(ny_out-1)+0.5);
            if (xx_out>nx_out-1)xx_out=nx_out-1;
            if (yy_out>ny_out-1)yy_out=ny_out-1;

			#if QSL_DIM==3
				k=zz_out+nz_out*(yy_out+ny_out*xx_out);
			#endif
			#if QSL_DIM==2
				k=(yy_out+ny_out*xx_out);
			#endif
			
			if (arr[k].qsl < q){//take the maximum qsl 
				arr[k].qsl = q;
			}
			if (arr[k].gap==-1) arr[k].gap=0;

        }
        fclose(file);
    }
    for (size_t i=0;i<nout;i++){
        cout << arr[i].qsl << "\t" << arr[i].gap;   
        cout<< "\n";
    }
    
}
            
