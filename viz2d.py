#
#	This file is part of QSL Squasher. 
#	Copyright (C) 2014, 2015, 2016  Svetlin Tassev
#							 Harvard-Smithsonian Center for Astrophysics
#							 Braintree High School
#	
#    QSL Squasher is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#   
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#   
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#   
#

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

########################################################################
########################################################################
# Uses the output from snapshot.cpp
########################################################################
########################################################################
arr=pd.read_table('grid2d.dat',header=None)

########################################################################
########################################################################
# Must match the corresponging values in options.hpp:
########################################################################
########################################################################
slice_center = [0,7.15]
slice_lx = 18.
slice_ly = 14.
CALCULATE='QSL'  # 'FIELD_LINE_LENGTH' or 'QSL'
########################################################################
########################################################################
# Must match the corresponging value in snapshot.cpp: 
########################################################################
########################################################################
nx_out=512*8
########################################################################
########################################################################
########################################################################




nx=nx_out
ny=len(arr)//nx
a=np.array(arr)[...,0].reshape((nx,ny)).T
xmin=slice_center[0]-slice_lx/2.
xmax=slice_center[0]+slice_lx/2.
ymin=slice_center[1]-slice_ly/2.
ymax=slice_center[1]+slice_ly/2.

if CALCULATE=='QSL':
    plt.imshow(a,vmin=-0.4,vmax=8,origin='lower',cmap=plt.get_cmap('jet'),extent=(xmin,xmax,ymin,ymax));plt.colorbar();
    #plt.savefig('grid2d.png', dpi=900)
    plt.show()
else:
    import scipy.ndimage.filters
    Gx = scipy.ndimage.filters.sobel(a,axis=0)
    Gy = scipy.ndimage.filters.sobel(a,axis=1)
    
    plt.imshow(np.sqrt(Gx**2+Gy**2),vmax=100,origin='lower',cmap=plt.get_cmap('gray'),extent=(xmin,xmax,ymin,ymax));plt.colorbar();
    plt.show()

