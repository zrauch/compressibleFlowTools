# compFlowTools
A library of compressible flow tools like area mach relation, mass flow function, isentropic relations, etc.
Originally developed during AAE 538 (Air-breathing propulsion) at Purdue University. Adapted over the duration of my graduate coursework including
ME 510 (Gas Dynamics), AAE 537 (Hypersonic propulsion), and AAE 539 (Rocket Propulsion).

# how to use
1. clone this repository to a local location on your working machine
  a. clone to the location where you wish to use these modules (your working directory for your .py script)
  b. (recommended) clone this to a general location and add that location to your sys.path in your .py script
 
If using option b, add the following lines to your .py script
import sys
\n sys.path
\n sys.path.append('~/path/to/repository/compFlowTools')

2. Once you have the library in a working location, import the modules you wish to use in your .py script
--> e.g "from compressible import *" or "import compressible as comp"

3. You can now start to use the compFlowTools function library as a Python module for quick access to important equation solvers for your compressible flow courses.
