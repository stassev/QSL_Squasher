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


.. _example-section:

Worked-out example
==================

The default options in the code generate a *very low-resolution* 3d cube 
of Q values with the sample dataset included with the code. Running the code with those 
options using the commands listed below, generates the following 3d view 
in ``ParaView``, which includes a magnetogram, volumetric rendering of 
Q values, as well as traced field lines:


.. figure::  3d_image.png
   :align:   center
   :width: 15cm
   :height: 11.25cm

The command lines below show a typical sequence of commands to generate 
a 3d data cube of Q values, and then visualize the result. The example 
uses the data files distributed with the code and was run on the *CPU* 
of a low-end laptop. The calculation of about half a million 
Q values, post-processing and rendering took less than 5 minutes with the 
default settings. Note that in the session below, qslSquasher was 
killed after two mesh refinements. ::

	$ time ./compile.sh 

	  real	0m20.256s
	  user	0m19.788s
	  sys	0m0.393s
	  
	$ time ./qslSquasher > raw.dat
	  1. pthread-Intel(R) Core(TM) i7-8650U CPU @ 1.90GHz (Portable Computing Language)
	  
	  Calculating current ...
	  ... done.
	  Analyzing transverse FL motions ...
	  ... done.
	  Reading successful.
	  Integration step set at 0.022446 Mm
	  Initialization successful.
	  Number of field lines to be integrated in this mesh refinement step: 262144
	  Beginning FORWARD integration along field lines ...
	  # of computed field lines = 0 out of 242172 in mesh refinement: 0
	     ... skipping lines ...
	  # of computed field lines = 242166 out of 242172 in mesh refinement: 0
	  Beginning BACKWARD integration along field lines ...
	  # of computed field lines = 0 out of 242172 in mesh refinement: 0
	     ... skipping lines ...
	  # of computed field lines = 242166 out of 242172 in mesh refinement: 0
	  Quantities calculated successfully.
	  Starting sort...
	  ... done sorting.
	  Starting sort...
	  ... done sorting.
	  Number of field lines to be integrated in this mesh refinement step: 265381
	  Beginning FORWARD integration along field lines ...
	  # of computed field lines = 0 out of 265223 in mesh refinement: 1
	     ... skipping lines ...
	  # of computed field lines = 265221 out of 265223 in mesh refinement: 1
	  Beginning BACKWARD integration along field lines ...
	  # of computed field lines = 0 out of 265223 in mesh refinement: 1
	     ... skipping lines ...
	  # of computed field lines = 265222 out of 265223 in mesh refinement: 1
	  Quantities calculated successfully.
	  Starting sort...
	  ... done sorting.
	  Starting sort...
	  ... done sorting.
	  Number of field lines to be integrated in this mesh refinement step: 557436
	  Beginning FORWARD integration along field lines ...
	  # of computed field lines = 0 out of 556914 in mesh refinement: 2
	     ... skipping lines ...
	  # of computed field lines = 3897603 out of 3897606 in mesh refinement: 5
	  Quantities calculated successfully.
	  ^C
	  
	  real	19m46.221s
	  user	119m21.692s
	  sys	2m16.436s
	$ time ./snapshot > grid3d.dat
	  
	  real	0m12.897s
	  user	0m12.611s
	  sys	0m0.260s
	$ time python viz3d.py

	  real	0m19.614s
	  user	0m18.337s
	  sys	0m1.218s
	$ paraview --state=viz3d_paraview.pvsm
	

The example above is for input in cartesian coordinates. It generates 
several output files:

* :file:`raw.dat` contains the raw output from :download:`qslSquasher.cpp`.

* :file:`grid3d.dat` is the result of the first post-processing step done by :download:`snapshot.cpp`.

* :file:`Global_Quantities.vtr` is a VTK file, containing the 
  rectilinear grid of global quantities in cartesian coordinates. 
  
* :file:`Local_Quantities.vtr` is a VTK file, containing the 
  rectilinear grid of local quantities in 
  cartesian coordinates. This file is generated from the input files 
  used by ``qslSquasher``.
  
The last two files are used by the included :download:`ParaView session file <viz3d_paraview.pvsm>` 
to generate the figure shown in the beginning of this section.
