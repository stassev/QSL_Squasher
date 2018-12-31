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


.. _options-section:

Options
========

The options for ``qslSquasher`` below are configurable in 
:download:`options.hpp`. This means that the code needs to be recompiled 
after each modification to the options below.


.. cpp:member:: std::string in_dir

   The directory name for the input files. See :ref:`Input <input-section>`.

.. cpp:member:: std::string in_filename

   The filename suffix for the input files. See :ref:`Input <input-section>`.


.. c:macro:: N(X|Y|Z)

   The size of the input arrays in the x/longitudinal (``NX``), 
   y/latitudinal(``NY``) and z/radial (``NZ``) direction.

.. c:macro:: GEOMETRY

   Pick one type of geometry for your input box. Possible values are: 
   ``SPHERICAL`` (default) or ``CARTESIAN``. When using spherical geometry,
   the poles, as well as the periodicity in longitude, are treated correctly.

.. c:macro:: GLOBAL_MODEL

   If defined (only for ``SPHERICAL`` geometry), the code assumes that 
   the input magnetic field covers the whole sun. Longitude samples should 
   start at >=0 degrees, and end at <360 degrees. Latitude samples should 
   start >-90 degrees, and end at <90 degrees. The code currently supports 
   only trilinear interpolation when this options is set. The poles, as well
   as the periodicity in longitude, are treated correctly.

.. c:macro:: PERIODIC_XY

   If set, then use periodic boundary conditions in X and Y.

.. c:macro:: SOLAR_RADIUS

   The solar radius in Mm (default: ``696.``). Needed only for 
   spherical geometry.
   
.. c:macro:: BOX_SIZE

   The typical size of the box in Mm. Used to set sane default values for 
   some of the other parameters of the code.

.. c:macro:: QSL_DIM
   
   Tells the code whether you want a 2d slice (``QSL_DIM=2``) of Q values, 
   or a 3d cube (``QSL_DIM=3``)?

If ``QSL_DIM=2``, then one needs to specify the parameters controlling 
the size, location and orientation of the slice one wants computed. 
Here is the set of relevant options that need to be set:
 
 .. c:macro:: SLICE_TYPE
 
	Specifies the type of slice. For cartesian geometry, the only 
	available option is ``CARTESIAN``. For spherical geometry, one can 
	pick a ``CARTESIAN`` or a ``SPHERICAL`` slice. When set to 
	``CARTESIAN``, the slice is an intersection the volume with a plane 
	of position, size and orientation specified by the options below. 
	When set to ``SPHERICAL`` (default), the slice is a curved 2d 
	surface at fixed radius.
 
 .. c:member:: double SLICE_NORMAL[]
 
	Vector normal to slice. Need not be normalized. Used only for cartesian slices.
 
 .. c:member:: double SLICE_UP[]
 
	``SLICE_UP`` gives the general "up" direction along the slice. We 
	take only the component of ``SLICE_UP`` that lies in the plane of 
	the slice to construct the ``y`` direction in the plane of the 
	slice. So, need not be orthonormal to ``SLICE_NORMAL``. Note that 
	the ``x`` direction in the plane of the slice is given by the cross 
	product ``SLICE_UP`` :math:`\times` ``SLICE_NORMAL``. So, be 
	careful with the overall sign of ``SLICE_UP``, or you may end up 
	with a flipped image. Used only for cartesian slices.
 
 .. c:member:: double SLICE_CENTER[]
 
	``SLICE_CENTER`` gives the coordinates of the center of the slice. 
	The coordinates are in units of (Mm, Mm, Mm) for cartesian 
	geometry, or in units of (degrees, degrees, Mm above the 
	photosphere) for spherical geometry. The slice will pass through 
	this point. 
 
 .. c:member:: double SLICE_L(X|Y)
 
	``SLICE_LX`` and ``SLICE_LY`` give the size of the slice in Mm for 
	cartesian slices. For spherical slices, the units are in degrees.
 
 .. c:macro:: ZMIN   
 
	``ZMIN`` forces field lines to be terminated at that height above 
	the photosphere/bottom of the box for spherical/cartesian 
	coordinates. This is useful for eliminating photospheric "noise".

If ``QSL_DIM=3``, then one needs to specify the size and location of 
the 3d cube for which the Q values are to be computed. Here is the set 
of relevant options that need to be set:

 .. c:macro:: (X|Y|Z)(MIN|MAX)

	These six parameters give the boundaries of the cube for the 3d Q 
	calculation. For cartesian geometry, all are in Mm. For spherical 
	geometry, the X and Y limits are set in degrees along the 
	longitudinal and latitudinal directions, respectively. In that 
	case, ``ZMIN`` and ``ZMAX`` are measured in Mm above the 
	photosphere. ``ZMIN`` also forces the calculation of the field 
	lines to terminate at that height above the solar 
	photosphere/bottom of the input box for spherical/cartesian 
	geometries. This is useful for eliminating photospheric "noise". 
	Note that apart from the ``ZMIN`` limit, the field lines are 
	followed to the boundaries of the data cube spanned by the *input* 
	files.

 .. c:macro:: z_sampler(z)
 
	A function specifying how to sample the 3d cube in the radial/z 
	direction for spherical/cartesian geometries. Its argument is 
	assumed normalized between ``0`` (corresponding to bottom index of 
	the cube) and ``1`` (corresponding to top index of the cube). Its 
	output must span the physical size of the box in Mm, i.e. it should 
	run between ``ZMIN`` and ``ZMAX``. 

.. c:macro:: CALCULATE

   The code calculates the local and global quantities associated with 
   the transverse eigenvalues of the gradient of the normalized 
   magnetic field when ``CALCULATE`` is set to 
   ``TRANSVERSE_EIGENVALUES`` (default). As part of the output, the 
   code will generate numerous .dat files containing the local 
   quantities associated with the transverse eigenvalues. Those are 
   post-processed by the python script and are output into the VTK 
   files. When set to ``QSL``, the code calculates the squashing factor 
   of the field lines passing through each sampled point. One has to go 
   through the same post-processing pipeline, irrespective of the 
   option set by ``CALCULATE``. 

.. c:macro:: RHO_Z

   When ``CALCULATE`` is set to ``TRANSVERSE_EIGENVALUES``, one can set
   several options for how the squeezing rate is calculated by setting
   ``RHO_Z`` to one of ``OPT1_LAMBDA``, ``OPT2_LAMBDA`` and ``SYMM_LAMBDA``.
   For the first two, see the equations for options 1 and 2 for
   :math:`\rho_{\mathcal{Z}}` in Section 2.2.3 of [Coiling]_. For the 
   third (symmetrized) option, see Section 4.5 of the same paper.
   

.. c:macro:: n(x|y|z)_init

   The size of the initial grid (before mesh refinement) for which the 
   Q values are to be computed. ``nz_init`` is not needed if 
   ``QSL_DIM=2``.


.. c:macro:: OpenCL_DEVICE_TYPE
   
   Tells VexCL whether to use the CPU when defined as 
   ``CL_DEVICE_TYPE_CPU`` (default), or the GPU when defined as 
   ``CL_DEVICE_TYPE_GPU``.


.. c:macro:: NGPU

   ``NGPU`` (default: ``0``) tells VexCL on which GPU you want to do 
   the computation. In case you want to specify the GPU in other ways, 
   consider changing the GPU filter specified by the following line in 
   :download:`qslSquasher.cpp`::
   
    vex::Context ctx(   vex::Filter::Type(OpenCL_DEVICE_TYPE) 
                     && vex::Filter::Position(NGPU) );



.. c:var:: const size_t  CHUNKSIZE
   
   The ``CHUNKSIZE`` sets how many Q value calculations are to be 
   dispatched to the GPU in one go. Set ``CHUNKSIZE`` too high and you'll 
   run out of GPU memory. Set it too low, and you'll find performance 
   being degraded. The proper value will depend mostly on your hardware 
   and on your choice for integration sheme, so experiment until you 
   find the sweet spot for your configuration. The default value 
   :math:`(2^{19})` is optimized for the ``EULER`` scheme on AMD FirePro 
   W8100, which has 8GB memory.


.. c:macro:: LENGTH_JUMP_REFINEMENT_THRESHOLD

   Specifies the threshold (default: 1Mm) for the change in field-line 
   length between two neighbouring points on the Hilbert curve. If that 
   threshold is exceeded, then the code makes a refinement by sampling 
   the point lying half-way on the Hilbert curve between those two 
   points.
   
.. c:macro:: MAX_REFINEMENTS

   Specifies the maximum number of refinements the code will make before
   exiting.


.. c:macro:: MAX_BATCHSIZE

   Specifies the maximum number of field lines to be integrated in each 
   refinement before the code exits.

.. c:macro:: LOCAL_Q

   By default, the Q value of a field line is obtained by calculating 
   the squashing factor between the two ends of a field line. An end of 
   a field line is considered the point where the field line intersects 
   the surface of the *input* b-field box, or where it hits a null.
   
   However, you can calculate a more localized value of Q by measuring 
   the squashing factor over a specified length (in Mm) up and down 
   each field line. To do that, uncomment the ``LOCAL_Q`` line in 
   :download:`options.hpp`. You'd need to specify the length over which 
   you want the local Q to be calculated. That is given by 
   ``INTEGRATION_RANGE`` in Mm. 

.. c:macro:: INTEGRATION_RANGE

   If ``LOCAL_Q`` is not defined, then the global Q values are 
   computed. In that case, field line integration is done in chunks 
   until the field line terminates at the box boundaries, or 
   :math:`\bm{B}` gets very close to zero (e.g. near nulls). The length 
   of each chunk is specified by ``INTEGRATION_RANGE``. After each 
   chunk, the field lines are checked for whether they have terminated. 
   If left undefined, a sane value for ``INTEGRATION_RANGE`` is picked 
   in :download:`qslSquasher.cpp`.

.. c:macro:: INTEGRATION_STEPS_PER_CELL

   ``INTEGRATION_STEPS_PER_CELL`` is used to calculate the step size 
   for the field line integrators. The step size is such that there are 
   roughly ``INTEGRATION_STEPS_PER_CELL`` steps in each cell in the 
   input grids. The resulting step size (printed to ``stderr``) is the 
   integration step for the Euler integration scheme. If left undefined, 
   sane defaults are set.
   

   
