.. ########################################################################
.. ########################################################################
.. #   This file is part of QSL Squasher. 
.. #   Copyright (C) 2014, 2015, 2016  Svetlin Tassev
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
   ``SPHERICAL`` (default) or ``CARTESIAN``. 

.. c:macro:: SOLAR_RADIUS

   The solar radius in Mm (default: ``696.``). Needed only for 
   spherical geometry.
   


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

   The code calculates the squashing factor values when ``CALCULATE`` 
   is set to ``QSL`` (default). When set to ``FIELD_LINE_LENGTH``, it 
   calculates the length of the field lines passing through each 
   sampled point. The code does not do refinements in the latter case, 
   as those are unnecessary for the field-line length map (as long as 
   the initial grid sampling is fine enough to resolve the connectivity 
   domains of interest). When calculating field-line lengths, the code 
   reuses the same infrastructure as when calculating the squashing 
   factor values. Thus, one has to go through the same post-processing 
   pipeline, irrespective of the option set by ``CALCULATE``.

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



.. c:macro:: INTERPOLATION_TYPE

   Pick one interpolation algorithm used for interpolating the B-field 
   values. Possible values are: ``TRILINEAR`` (default), ``TRIQUADRATIC``, 
   ``TRICUBIC``.
   
.. c:macro:: LENGTH_JUMP_REFINEMENT_THRESHOLD

   Specifies the threshold (default: 1Mm) for the change in field-line 
   length between two neighbouring points on the Hilbert curve. If that 
   threshold is exceeded, then the code makes a refinement by sampling 
   the point lying half-way on the Hilbert curve between those two 
   points.

.. c:macro:: INTEGRATION_SCHEME

   Specifies the integration scheme. One can set this to ``EULER`` 
   (default) for an explicit Euler scheme, or to ``ADAPTIVE`` for 
   adaptive stepping. The default adaptive stepper is Boost's 5-th 
   order `runge_kutta_cash_karp54 
   <http://www.boost.org/doc/libs/1_60_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html>`_. 
   You can always experiment with others by changing the corresponding 
   line in :download:`qslSquasher.cpp`.

.. c:macro:: eps_rel, eps_abs, DISPLACEMENT_WEIGHT

   Have an effect only when one uses the ``ADAPTIVE`` integration 
   scheme. The first two specify the relative and absolute error 
   (defaults: ``1e-2``) for the adaptive stepper. Those bounds are both 
   for the field line positions, as well as for the perturbations to the 
   field lines that are needed for the squashing factor calculation. 
   The ``DISPLACEMENT_WEIGHT`` (default: 10) boosts the weight of those 
   perturbations, since their errors will otherwise be swamped by the errors in 
   the positions.

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
   integration step for the Euler integration scheme, or is the initial 
   step for the adaptive stepper. If left undefined, sane defaults are 
   set in :download:`qslSquasher.cpp`.
   
.. c:macro:: MARK_OPEN_FIELD_LINES

   When ``MARK_OPEN_FIELD_LINES`` is defined (default), then the code 
   calculates `Q` values only for field lines which begin and end at 
   the bottom surface of the volume, corresponding to `z` or height 
   above the photosphere equal to `ZMIN` for cartesian or spherical 
   geometry, respectively. Open field lines are marked with the generic 
   value of `-1000`, which is used for any junk values encountered by the 
   code. If this keyword is left undefined, then the code calculates 
   the `Q` value for all points in the volume, irrespective of whether 
   they belong to open field lines or not. For open field lines, the 
   `Q` value is calculated between the two endpoints of the field 
   lines, independent of whether those occur at the bottom boundary or 
   not.
   
