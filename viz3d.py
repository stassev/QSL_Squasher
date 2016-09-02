#  
#  This file is part of QSL Squasher. 
#  Copyright (C) 2014, 2015, 2016  Svetlin Tassev
#  						 Harvard-Smithsonian Center for Astrophysics
#  						 Braintree High School
#  
#   QSL Squasher is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#  
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  

import numpy as np
from evtk.hl import gridToVTK
import pandas as pd
from numpy import exp

########################################################################
########################################################################
# Uses the output from snapshot.cpp
########################################################################
########################################################################
arr=pd.read_table('grid3d.dat',header=None)


########################################################################
########################################################################
# Must match corresponding definitions in options.hpp
########################################################################
########################################################################
solar_radius=696.

xmin   = -9.9    
xmax   =  9.9    
ymin   = -14.9   
ymax   =  14.9   
zmin   =  0.00001
zmax   =  14.3   
def z_sampler(z1):
    z=float(z1)
    z=z*(zmax-zmin)+zmin
    #z=(exp(((70.+z*300.)/370.)**3)-exp((70./370.)**3))/(exp(1.)-exp((70./370.)**3))*(zmax-zmin)+zmin
    return z

GEOMETRY="cartesian"
CALCULATE='QSL'  # 'FIELD_LINE_LENGTH' or 'QSL'

in_dir='./cartesian_demo/'
in_filename=''

########################################################################
########################################################################
# Must match corresponding definitions in snapshot.cpp
nx_out=128;
ny_out=128;
nz_out=128;

########################################################################
########################################################################


Lx      = xmax-xmin
Ly      = ymax-ymin

xx=((np.array(range(nx_out),dtype='float32'))/(nx_out-1.0)*Lx+xmin)
yy=((np.array(range(ny_out),dtype='float32'))/(ny_out-1.0)*Ly+ymin)

z_sampler=np.vectorize(z_sampler)
zz=z_sampler(np.array(range(nz_out),dtype='float32')/(nz_out-1))
log10q=(np.array(arr))[...,0].reshape((nx_out,ny_out,nz_out))

if (CALCULATE=='FIELD_LINE_LENGTH'):
    import scipy.ndimage.filters
    Gx = scipy.ndimage.filters.sobel(log10q,axis=0)
    Gy = scipy.ndimage.filters.sobel(log10q,axis=1)
    Gz = scipy.ndimage.filters.sobel(log10q,axis=2)
    log10q=np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)


if GEOMETRY=="cartesian":
	########################################################################
	# Export log10(Q) in cartesian coordinates:
	########################################################################
	

	#Need to .copy() the arrays before feeding to gridToVTK as otherwise you will get issues.
	gridToVTK("SquashingFactor_CartCoo",xx.copy(),yy.copy(),zz.copy(), pointData = {"log10(Q)" : (log10q.copy())} )
	
	
	########################################################################
	# Export B field to VTK file in cartesian coordinates:
	########################################################################
	
	cx=np.array(pd.read_table(in_dir+'xs0'+in_filename+'.dat',header=None))[...,0]
	cy=np.array(pd.read_table(in_dir+'ys0'+in_filename+'.dat',header=None))[...,0]
	cz=np.array(pd.read_table(in_dir+'zs0'+in_filename+'.dat',header=None))[...,0]
	
	nx=cx.size
	ny=cy.size
	nz=cz.size
	
	bx =   (np.array(pd.read_table(in_dir+'bx0'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	by =   (np.array(pd.read_table(in_dir+'by0'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	bz =   (np.array(pd.read_table(in_dir+'bz0'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	
	
	gridToVTK("MagneticField_CartCoo", cx.copy(), cy.copy(), cz.copy(), pointData = {"b" : (bx.copy(),by.copy(),bz.copy())})


if GEOMETRY=="spherical":
	
	# This is no longer a rectilinear grid in cartesian coordinates, so 
	# we need to specify the coordinates of each grid point.
	zzz=np.repeat(np.repeat([[zz]],ny_out,axis=1),nx_out,axis=0)
	xxx=np.repeat(np.repeat([[xx]],ny_out,axis=1),nz_out,axis=0).transpose((2,1,0))
	yyy=np.repeat(np.repeat([[yy]],nx_out,axis=1),nz_out,axis=0).transpose((1,2,0))
	cx=(zzz+solar_radius)*np.cos(xxx*np.pi/180.)*np.cos(yyy*np.pi/180.)
	cy=(zzz+solar_radius)*np.sin(xxx*np.pi/180.)*np.cos(yyy*np.pi/180.)
	cz=(zzz+solar_radius)*np.sin(yyy*np.pi/180.)
	
	mx=np.mean(cx)
	my=np.mean(cy)
	mz=np.mean(cz)
	#print "The center will be shifted to the origin for easier viz in ParaView."
	print "The center coordinates are: ",mx,my,mz
	
	#cx=cx-mx
	#cy=cy-my
	#cz=cz-mz
	# del xxx,yyy,zzz
	
	gridToVTK("SquashingFactor_SphCoo", cx.copy(), cy.copy(), cz.copy(), pointData = {"log10(Q)" : (log10q.copy())})


	########################################################################
	# Export B field to VTK file:
	########################################################################
	
	lons=np.array(pd.read_table(in_dir+'xs0'+in_filename+'.dat',header=None))[...,0]
	lats=np.array(pd.read_table(in_dir+'ys0'+in_filename+'.dat',header=None))[...,0]
	rads=np.array(pd.read_table(in_dir+'zs0'+in_filename+'.dat',header=None))[...,0]
	
	
	nx=lons.size
	ny=lats.size
	nz=rads.size
	
	#minus in b_theta is due to Griffiths' def
	b_phi   =   (np.array(pd.read_table(in_dir+'bx0'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	b_theta = - (np.array(pd.read_table(in_dir+'by0'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	b_r     =   (np.array(pd.read_table(in_dir+'bz0'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	
	r=np.repeat(np.repeat([[rads]],ny,axis=1),nx,axis=0)*solar_radius
	phi=np.repeat(np.repeat([[lons]],ny,axis=1),nz,axis=0).transpose((2,1,0))*np.pi/180.
	theta=np.pi/2.0-np.repeat(np.repeat([[lats]],nx,axis=1),nz,axis=0).transpose((1,2,0))*np.pi/180. # Griffiths' def
	cx=(r)*np.cos(phi)*np.sin(theta) #- mx
	cy=(r)*np.sin(phi)*np.sin(theta) #- my
	cz=(r)*np.cos(theta)             #- mz
	
	bx=b_r*np.sin(theta)*np.cos(phi) + b_theta*np.cos(theta)*np.cos(phi) + b_phi *(-np.sin(phi))
	by=b_r*np.sin(theta)*np.sin(phi) + b_theta*np.cos(theta)*np.sin(phi) + b_phi *(np.cos(phi))
	bz=b_r*np.cos(theta) + b_theta*(-np.sin(theta)) 
	
	
	gridToVTK("MagneticField_SphCoo", cx.copy(), cy.copy(), cz.copy(), pointData = {"b" : (bx.copy(),by.copy(),bz.copy())})
