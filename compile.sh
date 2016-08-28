#!/bin/bash

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

#########################
# Compile qslSquasher.cpp
#########################
clang -I/opt/intel/opencl-sdk/include qslSquasher.cpp -std=c++11 -lstdc++ -lm -I/usr/local/include/vexcl -lOpenCL -lboost_system -O3 -march=native -mcpu=native  -o qslSquasher
#clang  -I/opt/AMDAPP/SDK/include  qslSquasher.cpp -std=c++11 -lstdc++ -lm -I/usr/local/include/vexcl  -lOpenCL -lboost_system -O3 -o qslSquasher


######################
# Compile shapshot.cpp
######################
clang++  -std=c++11 snapshot.cpp -O3 -o snapshot  -march=native

##################################
# To build the documentation, run:
##################################
# sphinx-build -b html . _build/
