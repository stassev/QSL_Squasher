/*
 * This code calculates the Hilbert-Peano key and its inverse in 2D and 3D. 
 * The code is adopted from the locations identified below.
 * I'm licensing any changes to the original codes under GPLv3 or
 * (at your option) any later version of the GPL.
 */ 

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


/* 2D Hilbert code is written in C++ based on 
 * GeometricalPredicates.jl from https://gist.github.com/skariel/da85943803a6f57a52fd
 * Author: Ariel Keselman (skariel@gmail.com)
 * Original License: MIT 
 */ 

/* 3D Hilbert code is adopted from Gadget-2. Here is the original licence:
 * Copyright (c) 2005       Volker Springel
 *                          Max-Plank-Institute for Astrophysics
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <iostream>


#if QSL_DIM==3
	typedef  long long  peanokey_type;    /*!< defines the variable type used for Peano-Hilbert keys */
	
	
	//typedef long long uint64_t;
	
	static int quadrants[24][2][2][2] = {
	/* rotx=0, roty=0-3 */
	{{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
	{{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
	{{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
	{{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
	/* rotx=1, roty=0-3 */
	{{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
	{{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
	{{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
	{{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
	/* rotx=2, roty=0-3 */
	{{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
	{{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
	{{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
	{{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
	/* rotx=3, roty=0-3 */
	{{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
	{{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
	{{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
	{{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
	/* rotx=4, roty=0-3 */
	{{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
	{{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
	{{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
	{{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
	/* rotx=5, roty=0-3 */
	{{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
	{{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
	{{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
	{{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
	};
	
	
	static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
	12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
	};
	
	static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
	11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
	};
	
	static int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
	static int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };
	
	static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };
	
	static int flag_quadrants_inverse = 1;
	static char quadrants_inverse_x[24][8];
	static char quadrants_inverse_y[24][8];
	static char quadrants_inverse_z[24][8];
	
	

	uint64_t peanokey (uint64_t x, uint64_t y, uint64_t z)
	{
	int i, quad, bitx, bity, bitz;
	int mask, rotation, rotx, roty, sense;
	peanokey_type key;
	
	
	mask = 1 << (MAX_RES_BITS - 1);
	key = 0;
	rotation = 0;
	sense = 1;
	
	
	for(i = 0; i < MAX_RES_BITS; i++, mask >>= 1)
		{
		bitx = (x & mask) ? 1 : 0;
		bity = (y & mask) ? 1 : 0;
		bitz = (z & mask) ? 1 : 0;
	
		quad = quadrants[rotation][bitx][bity][bitz];
	
		key <<= 3;
		key += (sense == 1) ? (quad) : (7 - quad);
	
		rotx = rotx_table[quad];
		roty = roty_table[quad];
		sense *= sense_table[quad];
	
		while(rotx > 0)
		{
		rotation = rotxmap_table[rotation];
		rotx--;
		}
	
		while(roty > 0)
		{
		rotation = rotymap_table[rotation];
		roty--;
		}
		}
	
	return (uint64_t)key;
	}
	
	
	void point(uint64_t keyy, uint64_t *x, uint64_t *y, uint64_t *z)
	{
	peanokey_type key=keyy;
	long long i, keypart, bitx, bity, bitz, mask, quad, rotation, shift;
	char sense, rotx, roty;
	
	if(flag_quadrants_inverse)
		{
		flag_quadrants_inverse = 0ULL;
		for(rotation = 0ULL; rotation < 24ULL; rotation++)
			for(bitx = 0ULL; bitx < 2; bitx++)
			for(bity = 0ULL; bity < 2; bity++)
				for(bitz = 0ULL; bitz < 2; bitz++)
				{
					quad = quadrants[rotation][bitx][bity][bitz];
					quadrants_inverse_x[rotation][quad] = bitx;
					quadrants_inverse_y[rotation][quad] = bity;
					quadrants_inverse_z[rotation][quad] = bitz;
				}
		}
	
	shift = 3ULL * (MAX_RES_BITS - 1ULL);
	mask = 7ULL << shift;
	
	rotation = 0ULL;
	sense = 1ULL;
	
	*x = *y = *z = 0ULL;
	
	for(i = 0ULL; i < MAX_RES_BITS; i++, mask >>= 3ULL, shift -= 3ULL)
		{
		keypart = (key & mask) >> shift;
	
		quad = (sense == 1ULL) ? (keypart) : (7ULL - keypart);
	
		*x = (*x << 1ULL) + quadrants_inverse_x[rotation][quad];
		*y = (*y << 1ULL) + quadrants_inverse_y[rotation][quad];
		*z = (*z << 1ULL) + quadrants_inverse_z[rotation][quad];
	
		rotx = rotx_table[quad];
		roty = roty_table[quad];
		sense *= sense_table[quad];
	
		while(rotx > 0ULL)
			{
			rotation = rotxmap_table[rotation];
			rotx--;
			}
	
		while(roty > 0ULL)
			{
			rotation = rotymap_table[rotation];
			roty--;
			}
		}
	}
#endif


#if QSL_DIM==2
    uint64_t peanokey(uint64_t x, uint64_t y) {
		uint64_t n = 1 << MAX_RES_BITS, p=0, s = n >> 1;
        uint64_t rx, ry;
        while (s!=0) {
            rx = (x & s) > 0;
            ry = (y & s) > 0;
            p += s * s * ((3 * rx) ^ ry);
            if (ry == 0) {
				if (rx == 1) {
					x = n-1-x;
					y = n-1-y;
				}
				uint64_t tmp=x;
				x = y;
				y = tmp;
			}
            s>>=1;
        }
        return p;
    }
     
    void point(uint64_t peanokey, uint64_t *x, uint64_t *y) {
        uint64_t n = 1 << MAX_RES_BITS, s=1, rx, ry;
        *x=0, *y=0;
        while (s < n){
        
			rx = 1 & (peanokey >> 1);
			ry = 1 & (peanokey ^ rx);
			
			if (ry == 0){
				if (rx == 1){
					*x = s-1-*x;
					*y = s-1-*y;
				}
				uint64_t tmp  = *x;
				*x = *y;
				*y = tmp;
			}
			
			*x += s * rx;
			*y += s * ry;
					
			s <<= 1;
			peanokey >>= 2;
		}
        
    }
     
#endif


//int main1( int argc , char **argv ){
//
// 	uint64_t xx=5430;
//	uint64_t yy=662;
//	uint64_t zz=9843;
//	uint64_t pp= peanokey (xx,yy); //,zz);
//	printf("hello %lu %d %d %d \n",pp,xx,yy,zz);
//	xx=0;
//	yy=0;
//	zz=0;
//	point (pp,&xx,&yy) ; // ,&zz);
//	
//	printf("hello %lu %lu %d %d \n",pp,xx,yy,zz);
//	return 0;
//}
