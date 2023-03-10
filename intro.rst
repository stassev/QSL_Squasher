.. ########################################################################
.. ########################################################################
.. #   This file is part of QSL Squasher. 
.. #   Copyright (C) 2014-2019  Svetlin Tassev
.. #   						 Harvard-Smithsonian Center for Astrophysics
.. #   						 Braintree High School
.. #   
.. #    QSL Squasher is free software: you can redistribute it and/or modify
.. #    it under the terms of the GNU General Public License as published by
.. #    the Free Software Foundation, either version 3 of the License, or
.. #    (at your option) any later version.
.. #   
.. #    This program is distributed in the hope that it will be useful,
.. #    but WITHOUT ANY WARRANTY; without even the implied warranty of
.. #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. #    GNU General Public License for more details.
.. #   
.. #    You should have received a copy of the GNU General Public License
.. #    along with this program.  If not, see <http://www.gnu.org/licenses/>.
.. #   
.. ########################################################################
.. ########################################################################


Overview
============



QSL Squasher is an `OpenCL <https://www.khronos.org/opencl/>`_ code for 
calculating the squashing (Q) and squeezing (Z) factors as well as twist and 
coiling numbers of a vector field specified within 
a finite volume. Its description below focuses on its application to 
solar magnetic fields, but the code itself is completely general.

Note that this is the manual for the QSL Squasher -- Transverse Version 
(corresponding to QSL Squasher version >=2.0). This version of the code is 
restricted to Euler integration scheme and trilinear interpolation in 
calculating Q. For higher-order integration and interpolation schemes in the
calculation of Q, use the original version of QSL Squasher (version <2.0).

QSL Squasher is based on the following papers: [QSL3d]_ and  [Coiling]_. We kindly ask 
you [#f1]_ to acknowledge both papers and its authors in any program or 
publication in which you use QSL Squasher (or  QSL Squasher -- Transverse Version) 
or a derivative of it.


.. rubric:: Footnotes

.. [#f1] We cannot *require* you, however, as we want QSL Squasher to be 
   GPLv3+ compatible.

Compiling
---------


QSL Squasher requires `Boost <http://www.boost.org/>`_, `VexCL 
<https://github.com/ddemidov/vexcl>`_, an `OpenCL 
<https://www.khronos.org/opencl/>`_ implementation, as well as their 
respective dependencies. The visualization scripts require `Python 
<https://www.python.org/>`_ with `SciPy <https://www.scipy.org/>`_ and 
`PyEVTK <https://bitbucket.org/pauloh/pyevtk>`_. The resulting 3d data 
cubes are exported to `VTK <http://www.vtk.org/>`_ format, which can 
then be visualized using `Paraview <http://www.paraview.org/>`_, `VisIt 
<https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or `Mayavi 
<http://code.enthought.com/projects/mayavi/>`_ among many. 


The code has been mostly tested on an `Arch Linux 
<https://www.archlinux.org/>`_ server with an AMD FirePro W8100 GPU, 
and it has been successfully run on laptops with subpar hardware (the 
example included in this documentation was run on such a laptop). The 
code is memory hungry when performing large refinements in 3D, so we 
use a swap of 256GB on an SSD on the server. As a reference, the 
following relevant packages were installed on the server, which may or 
may not be required depending on your particular hardware 
configuration:


================== ==================
Package            Arch Linux Version
================== ==================
amd-adl-sdk        6.0-1
amdapp-aparapi     20130123-1
amdapp-codexl      1.6-7247
amdapp-sdk         2.9.1-1
amdapp-sdk-aparapi 2.9.1-1
amdapp-sdk-opencv  2.9.1-1
boost-compute-git  0.4-2
catalyst-firepro   14.502.1040-1
clang              3.6.2-2
intel-opencl-sdk   2014_R2-2
linux              4.1.5-1
pocl               0.11-1
vexcl-git          20150710-4
xorg-server        1.16.4-1 
================== ==================
	
A compile script, :download:`compile.sh`, is included with the source 
code which compiles the two main programs: the main calculation code, 
:download:`qslSquasher.cpp`, as well the post-processing code 
:download:`snapshot.cpp`. You need to edit the script by hand to make 
it consistent with your configuration, especially since different 
OpenCL implementations can co-exist on the same hardware on different 
paths.

As an example, for the AMD GPU on our dedicated server, 
:download:`qslSquasher.cpp` is compiled with:

.. code-block:: bash

   $ clang -I/opt/intel/opencl-sdk/include qslSquasher.cpp -std=c++11 \
   > -lstdc++ -lm -I/usr/local/include/vexcl -lOpenCL -lboost_system \
   > -O3 -march=native -mcpu=native  -o qslSquasher

If you want to test the code on your CPU, you would need the `POCL 
<http://www.portablecl.org>`_ OpenCL implementation to be installed on 
your computer. Then you need to define :c:macro:`OpenCL_DEVICE_TYPE 
<OpenCL_DEVICE_TYPE>` as ``CL_DEVICE_TYPE_CPU`` in 
:download:`options.hpp`. 


.. warning::

   Some of the options for the OpenCL kernels need to be specified at 
   compile time. Therefore, you need to recompile the code every time 
   you change the hard-coded options in :download:`options.hpp`. For 
   instance, by default, the code runs on the CPU, not on the GPU. To 
   run on the GPU, you need to set :c:macro:`OpenCL_DEVICE_TYPE 
   <OpenCL_DEVICE_TYPE>` to ``CL_DEVICE_TYPE_GPU``.


.. warning::

   Make sure you optimize the :c:data:`CHUNKSIZE` in 
   :download:`options.hpp` before using this code for production purposes. 
   If :c:data:`CHUNKSIZE` is set too high, you may run out 
   of GPU memory, and get curious error messages... 


.. _input-section:

Input
-----

The code is configured by adjusting the hard-coded values in 
:download:`options.hpp`. Those are described in detail :ref:`here 
<options-section>`.

The code takes as input 6 ASCII files with the following naming 
conventions::

    in_dir+'bx0'+in_filename+'.dat'
    in_dir+'by0'+in_filename+'.dat'
    in_dir+'bz0'+in_filename+'.dat'
    
    in_dir+'xs0'+in_filename+'.dat'
    in_dir+'ys0'+in_filename+'.dat'
    in_dir+'zs0'+in_filename+'.dat'


The files :file:`b(x|y|z)*.dat` contain the 3d arrays for 
the 3 components of the magnetic field as a flattened list of numbers. 
The units of the magnetic field can be arbitrary as only 
the tangent unit vectors are used to calculate the Q values. The 3d 
arrays are of dimensions (:c:macro:`NX, NY, NZ <N>`) and are read 
inside the following nested for loops::

	for (size_t k = 0; k < NZ; ++k)  
		for (size_t j= 0; j < NY; ++j) 
			for (size_t i = 0; i < NX; ++i)

So, take this ordering into account when writing inputs for this code.

The dimensions (:c:macro:`NX, NY, NZ <N>`) need to be set in 
:download:`options.hpp` at compile time for the OpenCL kernels, which 
means that you need to recompile the code for each new box. 


The code assumes that the magnetic field components are sampled on a 
rectilinear grid in either spherical or cartesian coordinates. The grid 
point coordinates are specified by the files :file:`xs*.dat`, 
:file:`ys*.dat`, :file:`zs*.dat`. Those samples should be in increasing
order.

Depending on whether :c:macro:`GEOMETRY <GEOMETRY>` is set to 
``CARTESIAN`` or ``SPHERICAL``, the magnetic field and the 
grid point coordinates are given as follows:

For the ``CARTESIAN`` setting, the files :file:`b(x|y|z)*.dat` contain 
the components of the magnetic field in the usual orthonormal 
:math:`\hat x,\ \hat y,\ \hat z` cartesian basis. For the ``SPHERICAL`` 
setting those files contain the magnetic field components in the 
orthonormal spherical basis :math:`\hat\phi,\ \hat\theta,\ \hat r` 
(i.e. longitude, latitude, radius). Therefore, in that case 
:file:`bx*.dat` contains the magnetic field at each grid point in the 
longitudinal direction, :file:`by*.dat` gives the latitudinal 
component, and :file:`bz*.dat` -- the radial. 

For the ``CARTESIAN`` setting, the grid coordinate files 
:file:`xs*.dat`, :file:`ys*.dat`, :file:`zs*.dat` contain :c:macro:`NX, 
NY, NZ <N>` numbers specifying the respective :math:`x,\ y,\ z` 
coordinates of the grid points in Mm. 

For the ``SPHERICAL`` setting, :file:`xs*.dat` contains :c:macro:`NX` 
numbers specifying the longitudes of the grid points in *degrees*, 
while the file :file:`ys*.dat` contains :c:macro:`NY` numbers 
specifying the latitudes of the grid points in *degrees*. The file 
:file:`zs*.dat` contains :c:macro:`NZ` numbers specifying the radial 
coordinates of the grid points in units of *solar radii*. In other 
words, the photosphere of the sun should be at :math:`r=1` in this file 
when using spherical coordinates.
  


Output
--------

.. _output-section:

Output from qslSquasher
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code calculates the squashing Q values for the input magnetic field 
on either a 2d slice or a 3d cube, depending on whether 
:c:macro:`QSL_DIM` is set to ``2`` or ``3``, respectively.

Progress and debugging information is output to ``stderr``, 
while the calculation results are output to ``stdout`` after the Q 
values are calculated for the initial grid, and then after each 
successive mesh refinement.

In the code, the slice/cube for which the Q values are calculated is 
indexed with a `Hilbert curve 
<https://en.wikipedia.org/wiki/Hilbert_curve>`_ that fills the region 
of interest. The output from the :download:`qslSquasher.cpp` code 
is printed to stdout. 
	
The output after the initial calculation on a grid and after each mesh 
refinement is sorted according to Hilbert key values. For multiple 
refinements, the output can easily reach more than a billion Q values 
sampled on an irregular grid. Thus, for convenience, we provide a 
series of post-processing routines, which allow for easier 
vizualization of the results. The post-processing is performed by 
:download:`snapshot.cpp` and the Python visualization scripts described 
below.

.. _snapshot-section:

First post-processing step with :file:`snapshot.cpp`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This code assumes that the output from ``qslSquasher`` is saved in the 
current directory as :file:`raw.dat`. Then, ``snapshot`` parses that 
file and returns to ``stdout`` a list of values for the quantities of interest
on a rectilinear grid spanning the 2d/3d region of interest, which was specified for
``qslSquasher``. The grid is of size ``nx_out``, ``ny_out`` 
(and ``nz_out`` when working with a 3d cube), which are specified at 
the top of :download:`snapshot.cpp`  at compile time. 

The output is a column of values printed by the following nested for-loops::

	 for (size_t i = 0; i < nx_out; ++i)
		for (size_t j= 0; j < ny_out; ++j) 
			for (size_t k = 0; k < nz_out; ++k) # for 3d cube

The output is parsed with the Python script described in the next 
section. See that section for a description of the output.


.. _script-section:

Second post-processing step with Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The script :file:`viz3d.py` shows examples of 
post-processing the 3d output from ``snapshot``. 

The 3d post-processing script, :file:`viz3d.py`, outputs two VTK files 
containing the global and local quantities in the 3d 
cube sampled by qslSquasher. Those VTK files can then be visualized by 
ParaView as in the :ref:`example <example-section>` included with this 
documentation.

The quantities saved in the VTK files by the Python script are as follows:

"Type": 0 for transverse saddle flow; 1 for node; 2 for LH spiral; 3 for RH spiral; 
4 for LH center; 5 for RH center; 6 for improper node (one repeated eigenvector);
7 for star node (two distinct eigenvectors). 

"rho_Z" = :math:`\rho_{\mathcal{Z}}`, the squeezing rate (units of 1/Mm)

"omega_c" = :math:`\omega_c`, the coiling rate (units of 1/Mm)

"Trace" = sum of transverse eigenvalues of the gradient of the normalized magnetic field.

"J" = current (defined as :math:`\nabla \times \vec{B}`; has units of [B]/Mm.)   

"Alpha" = generalized force-free parameter (units of 1/Mm)

"Alpha_im" = generalized force-free parameter, which is set to zero in regions with real tranverse eigenvalues (units of 1/Mm)

"FLL" = field-line length (in Mm)

"open" = 1 for open field lines; 0 for closed field lines; otherwise, for missing data

"N_c" = coiling number

"log10(Z)" = logarithm of the squeeze factor

"log10(Q)" = logarithm of the squashing factor

"N_t" = standard twist number

"N_t_im" = twist number after dropping saddle-flow contributions (see 
Section 5.3 of [Coiling]_)

"FLEDGE" = FLEDGE map

"open_dilat" = the "open" array with removed boundary pixels (between open and closed field-line regions); useful for 
filtering only closed field lines of the FLEDGE map, without including the large-FLEDGE boundary pixels.
