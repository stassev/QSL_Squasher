**Author:** Svetlin Tassev (Harvard-Smithsonian Center for Astrophysics, Braintree High School)

**Initial public release date:** Aug 28, 2016

QSL Squasher is an OpenCL code for calculating squashing factors of vector fields.
Version >=2.0 (aka QSL Squasher -- Transverse Version) can also calculate local and global quantities associated with the transverse eigenvalues of the gradient of the normalized magnetic field. Those include the squeeze factor and squeezing rate; the coiling number and coiling rate; the twist factor and twist factor constrained to non-saddle-flow regions; etc.

For efficient calculation of the squashing factor, use the original version of QSL Squasher (version <2.0). Please, use QSL Squasher -- Transverse Version (version >=2.0) for all calculations involving the squashing factor (when speed is not of essence) as well as the quantities associated with the transverse eigenvalues of the gradient of the normalized magnetic field.

QSL Squasher (<2.0) is based on the following paper:

* QSL Squasher: A Fast Quasi-Separatrix Layer Map Calculator, S. Tassev and A. Savcheva, [arXiv:1609.00724](https://arxiv.org/abs/1609.00724), The Astrophysical Journal, 840, 89 (2017).

QSL Squasher -- Transverse Version (version >=2.0) is based on the paper above as well as on:

* Coiling and Squeezing: Properties of the Local Transverse Deviations of Magnetic Field Lines, S. Tassev and A. Savcheva (2019), [arXiv:1901.00865](https://arxiv.org/abs/1901.00865).

If you use any version of QSL Squasher or a derivative of it for scientific work, we 
kindly ask you to reference the papers above.

* QSL Squasher is free and open-source software, distributed under the GPLv3+ license.

* To build the code and learn how to run it, read the manual [here](https://bitbucket.org/tassev/qsl_squasher/downloads/QSLSquasher.pdf) and [here](https://bitbucket.org/tassev/qsl_squasher/downloads/QSLSquasherTrans.pdf). Scripts that compile the code and its dependencies under CentOS and Debian are available from the author upon request.

* The example input files can be downloaded [here](https://bitbucket.org/tassev/qsl_squasher/downloads/cartesian_demo.tar.gz).

**Revision History:**

ver. 2.0 (Jan 1, 2019): QSL Squasher -- Transverse Version. Major overhaul of the code focusing on calculating local and global quantities associated with the transverse eigenvalues of the gradient of the normalized magnetic field. Those include: the squeeze factor as well as the coiling and twist numbers.

ver. 1.3 (Jan 22, 2018): 
	(1) Add support for periodic BC in X and Y for Cartesian coordinates. To enable, set PERIODIC_XY in options.hpp. 
	(2) Overhaul of snapshot.cpp. Gaps are interpolated with python instead. This update improves memory consumption and fixes certain artefacts due to the old Hilbert-curve based interpolation of the data. Thanks to Roger Scott (University of Dundee) for reporting those artefacts.

ver. 1.2 (Jan 5, 2018): Added terms neglected in eq.6 of paper. Terms have no effect on identification of QSL locations. Only Q values in spherical coordinates are affected. The extra contributions are suppressed by the ratio (\delta r/R), where \delta r is the scale over which B varies and R is the radius of the Sun. 

ver. 1.1 (June 27, 2017): Added support for global models (covering the Sun in longitude and latitude). Treat the poles as well as the periodicity in longitude correctly.
