----------------------------------------------------------------------
                                                                      
   A Library to Exchange Boundary Conditions via Sockets in OpenFOAM  
                                                                      
----------------------------------------------------------------------
 Blame : Dominic Chandar
 dominic.chandar@gmail.com




Currently Tested on OpenFOAM v2.3.0

Contents are as follows:

1.  OpenFOAM/db/socketStream: Intended to be in the "stream" format of C++, the current 
    version used standard send/recv for socket communication. Send and receive is 
    performed via special boundary conditions described below. Future releases will have
    the required operators overloaded for easy send and receive.

2.  finiteVolume/..../socketFvPatchField{V}: Boundary conditions exchange via sockets.
    Each patch initiates a socket connection. Port numbers change according to the variable.


Features and Drawbacks
-----------------------

*  Can run 2 OpenFOAM solvers on different machines (ONLY) in parallel - The serial version needs
   a small fix.

*  Currently runs the exact same solver on Server and Client. Users can however modify the code easily
   to run different solvers on the Server and Client respectively. The server and client communicate when
   correctBoundaryConditions() are called.

*  The server and the client both interpolate data using the pointToPointPlanarInterpolation class. However
   the number of boundary faces must match - A minor fix to make arbitrary number of faces.

*  The test case in tutorials/incompressible/simpleFoam demonstrates running simpleFoam for the flow past an 
   airfoil where a "near" body region uses a turbulence model whereas the "far" region does not.


