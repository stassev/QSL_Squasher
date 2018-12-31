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
  float fll;
  int ofl;
  #if CALCULATE==QSL
	float qsl;
  #else
	float a;
	float c;
	float b;
	float d;
  #endif
} ;

typedef std::vector< out_struct > out_type;


int main( void )
{
    using namespace std;
    uint64_t nin=0;
    out_type arr(nout);
    FILE *file;
	for (size_t kk=0;kk<nout;kk++){
		#if CALCULATE==QSL
			arr[kk].qsl=-1000;
		#else
		arr[kk].a=0;
		arr[kk].b=0;
		arr[kk].c=0;
		arr[kk].d=0;
		#endif
		arr[kk].fll=0;
		arr[kk].ofl=-1;
	}
    //cerr<<"ok0a\n";
    {
        float q,m,bx,by,bz,h,tmp,a,v,t,fll;
        uint64_t d,k;
        int ofl;
        file = fopen((raw_file).c_str(), "r");
		float b,c,dd;
        uint64_t xx,yy,zz;
        uint64_t xx_out,yy_out,zz_out;

        size_t i=0;
        
        while (!feof(file)) {
			nin+=1;
			#if CALCULATE!=QSL
				int check=fscanf(file, "%lu%g%d%g%g%g%g", &d, &fll, &ofl, &a, &b, &c, &dd );
				if ((check!=7)&&(check>0)){
			#else
			int check=fscanf(file, "%lu%g%g%d", &d, &q, &fll, &ofl);
				if ((check!=4)&&(check>0)){
			#endif
					std::cerr<<"***Problem with raw file at line "<<nin<<'\n';
					exit(0);
			    }
         #if CALCULATE==QSL
            if (q>1e-6) q=log10(q);       
            else q=-1000;                     
            
            if ((q>100) || (q<-6) || (q!=q)){// Check for junk values.
                q=-1000;
            }
            double qsl=q;
         #else
			if (dd!=dd) dd=0;
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
			
				#if CALCULATE==QSL
					if (arr[k].qsl < qsl){//take the maximum Q
						arr[k].qsl = qsl;
						arr[k].fll=fll;
					}
                #else
				if (
				    ((fabs(arr[k].b) == 0) && (fll>0) &&  (fabs(arr[k].a) < fabs(a))) 
				    || (fabs(arr[k].b) < fabs(b))
				    || ((fabs(arr[k].a) == 0) && (fabs(arr[k].b) == 0) && (fll>0))
				    ) {//taking maximum Z/Nc
					arr[k].d=dd;
					arr[k].c=c;
					arr[k].b=b;
					arr[k].a=a;
  				    arr[k].fll=fll;
					}
                #endif

				if (arr[k].ofl==-1) arr[k].ofl=0;
				if (ofl==1) arr[k].ofl=1;
        }
        fclose(file);
    }
    for (size_t i=0;i<nout;i++){
		#if CALCULATE==QSL
        cout << arr[i].fll << "\t" << arr[i].ofl<< "\t" << arr[i].qsl ;   
        #else
        cout << arr[i].fll << "\t" << arr[i].ofl<< "\t" << arr[i].a<< "\t" << arr[i].b<< "\t" << arr[i].c<< "\t" << arr[i].d ;   
        #endif
        cout<< "\n";
    }
    
}
            
