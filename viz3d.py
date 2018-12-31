#  
#  This file is part of QSL Squasher. 
#  Copyright (C) 2014-2019  Svetlin Tassev
#  						    Harvard-Smithsonian Center for Astrophysics
#  						    Braintree High School
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
from pyevtk.hl import gridToVTK
#from evtk.hl import gridToVTK # Choose this or the top line, depending on evtk version.
import pandas as pd
from numpy import exp

def read_crappy_unformatted_idl_file(filename):
	with open(filename,'r') as f:
		#next(f) # skip first row
		df = pd.DataFrame(l.rstrip().split() for l in f)
	dd=df.values.flatten()
	dd=dd[dd != np.array(None)]
	return np.array(dd).astype(np.float32)

#import numpy as np
from scipy import ndimage as nd

def fill_gaps(data, ofl): 
    mask=ofl*0
    mask[np.where(ofl<-0.5)]=True
    mask[np.where(np.isinf(data))]=True
    mask[np.where(np.fabs(data)>1.e5)]=True
    mask[np.where(np.isnan(data))]=True
    ind = nd.distance_transform_edt(mask, 
                                    return_distances=False, 
                                    return_indices=True)
    return data[tuple(ind)]

def erode_closed_field_line_region(a,width=2):
	struct = nd.generate_binary_structure(3, 3)
	return nd.binary_dilation(a, structure=struct,iterations=width).astype(a.dtype)
########################################################################
########################################################################
# Uses the output from snapshot.cpp
########################################################################
########################################################################
arr=pd.read_table('grid3d.dat',header=None,dtype=np.float32)

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

GEOMETRY="cartesian" #"spherical" or "cartesian"
CALCULATE='TRANSVERSE_EIGENVALUES'  # 'TRANSVERSE_EIGENVALUES' or 'QSL'

in_dir_base='./'
in_dir=in_dir_base+'cartesian_demo/'
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

if (CALCULATE=='QSL'):
	log10q=(np.array(arr))[...,2].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	ofl=(np.array(arr))[...,1].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	log10q[np.where(np.isinf(log10q))]=100
	ofl[np.where(log10q<0)]=-1
	log10q=fill_gaps(log10q,ofl)
	fll=(np.array(arr))[...,0].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	ofl=(np.array(arr))[...,1].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	ofl[np.where(fll<0)]=-1
	fll=fill_gaps(fll,ofl)
	ofl=fill_gaps(ofl,ofl.copy())
	ofl_dilat=erode_closed_field_line_region(np.int32(ofl+0.3)).astype(np.float32)
else:
	fll 		=(np.array(arr))[...,0].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	INTimLambda =(np.array(arr))[...,2].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	INTdLambda  =(np.array(arr))[...,3].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	INTalpha    =(np.array(arr))[...,4].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	INTalphaIm  =(np.array(arr))[...,5].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	gaps 		=(np.array(arr))[...,1].reshape((nx_out,ny_out,nz_out)).astype(np.float32)
	ofl=gaps.copy()
	ofl[np.where(fll<0)]=-1
	fll=fill_gaps(fll,ofl)
	ofl=fill_gaps(ofl,ofl.copy())
	ofl_dilat=erode_closed_field_line_region(np.int32(ofl+0.3)).astype(np.float32)
	
	INTimLambda =fill_gaps(INTimLambda ,gaps)
	INTalpha    =fill_gaps(INTalpha    ,gaps)
	INTalphaIm  =fill_gaps(INTalphaIm  ,gaps)
	
	INTimLambda/=(2.*np.pi)
	INTalpha/=(4.*np.pi)
	INTalphaIm/=(4.*np.pi)
	INTdLambda*=0.43429448 #log10(exp(1))
	INTdLambda+=0.30103 # log10(2) set ic Z=2
	INTdLambda[np.where(np.isinf(INTdLambda))]=100
	INTdLambda  =fill_gaps(INTdLambda  ,gaps)

####
####
####

import scipy.ndimage.filters
Gx = scipy.ndimage.filters.sobel(fll,axis=0)
Gy = scipy.ndimage.filters.sobel(fll,axis=1)
Gz = scipy.ndimage.filters.sobel(fll,axis=2)
fledge=np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz).astype(np.float32)


if GEOMETRY=="cartesian":
	########################################################################
	# Export log10(Q) in cartesian coordinates:
	########################################################################
	
	if (CALCULATE=='QSL'):
		gridToVTK(in_dir_base+"Global_QuantitiesQ",xx.astype(np.float32).copy(),yy.astype(np.float32).copy(),zz.astype(np.float32).copy(), pointData = 
			{"FLL" : (fll.astype(np.float32).copy()),
			 "open" : (ofl.astype(np.float32).copy()),
			 "log10(Q)" : (log10q.astype(np.float32).copy()),
			 "FLEDGE" : fledge.astype(np.float32).copy(),
			 "open_dilat" : ofl_dilat.copy()
			 }) 
		stop
	else:
		 gridToVTK(in_dir_base+"Global_Quantities",xx.astype(np.float32).copy(),yy.astype(np.float32).copy(),zz.astype(np.float32).copy(), pointData = 
			{"FLL" : (fll.astype(np.float32).copy()),
			 "open" : (ofl.astype(np.float32).copy()),
			 "N_c" : INTimLambda.astype(np.float32).copy(),
			 "log10(Z)" : INTdLambda.astype(np.float32).copy(),
			 "N_t" : INTalpha.astype(np.float32).copy(),
			 "N_t_im" : INTalphaIm.astype(np.float32).copy(),
			 "FLEDGE" : fledge.astype(np.float32).copy(),
			 "open_dilat" : ofl_dilat.copy()
			 }) 
	#stop
	
	########################################################################
	# Export B field to VTK file in cartesian coordinates:
	########################################################################
	
	cx=np.array(pd.read_table(in_dir+'xs0'+in_filename+'.dat',header=None))[...,0].astype(np.float32)
	cy=np.array(pd.read_table(in_dir+'ys0'+in_filename+'.dat',header=None))[...,0].astype(np.float32)
	cz=np.array(pd.read_table(in_dir+'zs0'+in_filename+'.dat',header=None))[...,0].astype(np.float32)
	
	nx=cx.size
	ny=cy.size
	nz=cz.size
	
	bx =   (read_crappy_unformatted_idl_file(in_dir+'bx0'+in_filename+'.dat').reshape((nz,ny,nx)) ).transpose((2,1,0)).astype(np.float32)
	by =   (read_crappy_unformatted_idl_file(in_dir+'by0'+in_filename+'.dat').reshape((nz,ny,nx)) ).transpose((2,1,0)).astype(np.float32)
	bz =   (read_crappy_unformatted_idl_file(in_dir+'bz0'+in_filename+'.dat').reshape((nz,ny,nx)) ).transpose((2,1,0)).astype(np.float32)

	nx-=1
	ny-=1
	nz-=1
	
	dLambda 	=   (np.array(pd.read_table(in_dir_base+'ReDeltaLambda.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Type 		=   (np.array(pd.read_table(in_dir_base+'ODE_type.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Alpha_im 	=   (np.array(pd.read_table(in_dir_base+'Alpha_im.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	ImLambda 	=   (np.array(pd.read_table(in_dir_base+'ImLambda.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Trace 		=   (np.array(pd.read_table(in_dir_base+'Trace.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Alpha     	=   (np.array(pd.read_table(in_dir_base+'Alpha.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	
	jx =   (np.array(pd.read_table(in_dir_base+'Jx'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0)).astype(np.float32)
	jy =   (np.array(pd.read_table(in_dir_base+'Jy'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0)).astype(np.float32)
	jz =   (np.array(pd.read_table(in_dir_base+'Jz'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0)).astype(np.float32)
	
	gridToVTK(in_dir_base+"Local_Quantities", cx.astype(np.float32).copy(), cy.astype(np.float32).copy(), cz.astype(np.float32).copy(), 
		cellData = 		
			{"Type" 	 : Type.astype(np.float32).copy(), 
			 "rho_Z"   : dLambda.astype(np.float32).copy(), 
			 "omega_c"  : ImLambda.astype(np.float32).copy(), 
			 "Trace" 	 : Trace.astype(np.float32).copy(), 
			 "J"         : (jx.astype(np.float32).copy(),jy.astype(np.float32).copy(),jz.astype(np.float32).copy()),
			 "Alpha_im"	 : Alpha_im.astype(np.float32).copy(),
			 "Alpha"     : Alpha.astype(np.float32).copy()
			 },
		pointData = {"b" : (bx.copy(),by.copy(),bz.copy())})
	
if GEOMETRY=="spherical":
	
	# This is no longer a rectilinear grid in cartesian coordinates, so 
	# we need to specify the coordinates of each grid point.
	zzz=np.repeat(np.repeat([[zz]],ny_out,axis=1),nx_out,axis=0).astype(np.float32)
	xxx=np.repeat(np.repeat([[xx]],ny_out,axis=1),nz_out,axis=0).transpose((2,1,0)).astype(np.float32)
	yyy=np.repeat(np.repeat([[yy]],nx_out,axis=1),nz_out,axis=0).transpose((1,2,0)).astype(np.float32)
	cx=(zzz+solar_radius)*np.cos(xxx*np.pi/180.)*np.cos(yyy*np.pi/180.)
	cy=(zzz+solar_radius)*np.sin(xxx*np.pi/180.)*np.cos(yyy*np.pi/180.)
	cz=(zzz+solar_radius)*np.sin(yyy*np.pi/180.)
	
	mx=np.mean(cx)
	my=np.mean(cy)
	mz=np.mean(cz)
	
	print("The center coordinates are: ",mx,my,mz)
	if (CALCULATE=='QSL'):
		gridToVTK(in_dir_base+"Global_QuantitiesQ",cx.astype(np.float32).copy(),cy.astype(np.float32).copy(),cz.astype(np.float32).copy(), pointData = 
			{"FLL" : (fll.astype(np.float32).copy()),
			 "open" : (ofl.astype(np.float32).copy()),
			 "log10(Q)" : (log10q.astype(np.float32).copy()),
			 "FLEDGE" : fledge.astype(np.float32).copy(),
			 "open_dilat" : ofl_dilat.copy(),	 
			 })  
		stop
	else:
		gridToVTK(in_dir_base+"Global_Quantities",cx.astype(np.float32).copy(),cy.astype(np.float32).copy(),cz.astype(np.float32).copy(), pointData = 
			{"FLL" : (fll.astype(np.float32).copy()),
			 "open" : (ofl.astype(np.float32).copy()),
			 "N_c" : INTimLambda.astype(np.float32).copy(),
			 "log10(Z)" : INTdLambda.astype(np.float32).copy(),
			 "N_t" : INTalpha.astype(np.float32).copy(),
			 "N_t_im" : INTalphaIm.astype(np.float32).copy(),
			 "FLEDGE" : fledge.astype(np.float32).copy(),
			 "open_dilat" : ofl_dilat.copy(),	 
			 })  
			 
   	#stop
	########################################################################
	# Export Local Quantities to VTK file:
	########################################################################
	
	lons=np.array(pd.read_table(in_dir+'xs0'+in_filename+'.dat',header=None))[...,0].astype(np.float32)
	lats=np.array(pd.read_table(in_dir+'ys0'+in_filename+'.dat',header=None))[...,0].astype(np.float32)
	rads=np.array(pd.read_table(in_dir+'zs0'+in_filename+'.dat',header=None))[...,0].astype(np.float32)
	
	
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
	
	bx=(b_r*np.sin(theta)*np.cos(phi) + b_theta*np.cos(theta)*np.cos(phi) + b_phi *(-np.sin(phi))).astype(np.float32)
	by=(b_r*np.sin(theta)*np.sin(phi) + b_theta*np.cos(theta)*np.sin(phi) + b_phi *(np.cos(phi))).astype(np.float32)
	bz=(b_r*np.cos(theta) + b_theta*(-np.sin(theta)) ).astype(np.float32)

	lons=((lons[0:nx-1]+lons[1:nx])/2.).copy()
	lats=((lats[0:ny-1]+lats[1:ny])/2.).copy()
	rads=((rads[0:nz-1]+rads[1:nz])/2.).copy()
		
	nx-=1
	ny-=1
	nz-=1
	
	r=np.repeat(np.repeat([[rads]],ny,axis=1),nx,axis=0)*solar_radius
	phi=np.repeat(np.repeat([[lons]],ny,axis=1),nz,axis=0).transpose((2,1,0))*np.pi/180.
	theta=np.pi/2.0-np.repeat(np.repeat([[lats]],nx,axis=1),nz,axis=0).transpose((1,2,0))*np.pi/180. # Griffiths' def
		
	dLambda 	=   (np.array(pd.read_table(in_dir_base+'ReDeltaLambda.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Type 		=   (np.array(pd.read_table(in_dir_base+'ODE_type.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Alpha_im 	=   (np.array(pd.read_table(in_dir_base+'Alpha_im.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	ImLambda 	=   (np.array(pd.read_table(in_dir_base+'ImLambda.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Trace 		=   (np.array(pd.read_table(in_dir_base+'Trace.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	Alpha     	=   (np.array(pd.read_table(in_dir_base+'Alpha.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	
	j_phi   =   (np.array(pd.read_table(in_dir_base+'Jx'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	j_theta = - (np.array(pd.read_table(in_dir_base+'Jy'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	j_r     =   (np.array(pd.read_table(in_dir_base+'Jz'+in_filename+'.dat',header=None))[...,0].reshape((nz,ny,nx)) ).transpose((2,1,0))
	
	jx=j_r*np.sin(theta)*np.cos(phi) + j_theta*np.cos(theta)*np.cos(phi) + j_phi *(-np.sin(phi))
	jy=j_r*np.sin(theta)*np.sin(phi) + j_theta*np.cos(theta)*np.sin(phi) + j_phi *(np.cos(phi))
	jz=j_r*np.cos(theta) + j_theta*(-np.sin(theta)) 
	
	gridToVTK(in_dir_base+"Local_Quantities", cx.astype(np.float32).copy(), cy.astype(np.float32).copy(), cz.astype(np.float32).copy(), cellData = 
			{"Type" 	 : Type.astype(np.float32).copy(), 
			 "rho_Z"   : dLambda.astype(np.float32).copy(), 
			 "omega_c"  : ImLambda.astype(np.float32).copy(), 
			 "Trace" 	 : Trace.astype(np.float32).copy(), 
			 "J"         : (jx.astype(np.float32).copy(),jy.astype(np.float32).copy(),jz.astype(np.float32).copy()),
			 "Alpha_im"	 : Alpha_im.astype(np.float32).copy(),
			 "Alpha"     : Alpha.astype(np.float32).copy()
			 },
			 pointData = {"b" : (bx.astype(np.float32).copy(),by.astype(np.float32).copy(),bz.astype(np.float32).copy()), "b_sph" : (b_r.astype(np.float32).copy(),b_theta.astype(np.float32).copy(),b_phi.astype(np.float32).copy())})
		
