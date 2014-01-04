PROTEUS -- PRotean Open-source magneTohydrodynamic Equations Unified Solver

Proteus : An early sea-god said to be capable of both telling the future as well as shape-shifting

Protean : Tending or able to change frequently or easily.

Proteus is an ambitious project to unify many of the numerical methods common to solving computational hydrodynamic and magnetohydrodynamic problems under a unifying framework. This is a work-in-progress, but at this juncture it is assumed if you are reading this file that you are really interested in using IMHD, which solved the incompressible equations of MHD, and is the progenitor from which Proteus is spawned.  Instructions below will be for using IMHD, and will most likely be invalid unless you check out the appropriate tag from the github repository.

0.  IMHD Overview
1.  Compilation
2.  Running
3.  Numerical Details
3a. -- Psuedo-Spectral algorithms
3b. -- Solenoidal Condition
3c. -- Misc Features

===============================================================================
0. IMHD Overview
-------------------------------------------------------------------------------

IMHD, which stands for Incompressible Magnetohydrodynamics, is one of the major codes I have used for my thesis work.  The name is actually a misnomer as the code actually solves a more versatile set of equations, but the code has primarily been used in the incompressible regime, where convection is disallowed, and so the name has stuck. In general, the code solves the following equations:

∇ · u = 0,      ∇ · B = 0 
∂_t u + ∇·(u u - B B) = -∇P + Ra Pr T + B^2 / β + Pr ∇^2 u + F_u 
∂_t u + ∇·(u T) = u·k̂ + ∇^2 T 
∂_t B = ∇ x (u x B) + Pr / Pm ∇^2 B + F_B

These are a fairly standard set of Boussinesq equations with a few minor additions.  We have the velocity field u and magnetic field B, the pressure P and the temperature field T.  Most of these terms should be familiar to people that work on MHD problems, but notable additions include F_u and F_B -- arbitrary forcing terms that have proven useful to drive various problems -- and an artificial magnetic buoyancy term B^2 / β.

Each term in these equations can be independently enabled or disabled to create a very customization set of equations. 

===============================================================================
1. Compilation
-------------------------------------------------------------------------------

There should be three directories present -- Run, Compile and src. The first and the second contain sample execution scripts and Makefiles, which you will need to configure yourself for usage on your local cluster, while the third contains the actual source files.

It generally does not take much effort to modify a makefile for a new computer. The only requirements this code has is that the local machine must have fftw3 installed (any version 3.x will probably work) and a mpi wrappers for c and c++ compilers. If you are compiling FFTW manually, it is best to compile both the single and double precision libraries.  The name of the c and c++ wrappers should be stored in the Makefile's cc and CC variables respectively. Any compiler options required to find and link to the FFTW3 installation should be contained in LIBS. If the FFTW3 header files are not on the automatic search path then they should be added to the CCFLAGS line as well.

The only two aspects of the code that should have to be altered at compile time, the first of which is the logging level. This can be changed by editing src/include/Log.h. Here one has the choice of setting the global code default to either Trace, Debug, Info, Warning, or Error level logging. By default this is set to Info by including LogInfo.h. Log.h can also be included anywhere in the code to locally change the logging level from that point to the end of file. This serves as a rough mechanism to allow detailed outputs in a controlled section of code only.

The other compile-time alteration possible is the level of precision the code operates in.  single precision allows the code to run faster, especially on network-bound systems where MPI will now only have to shuffle half as many bits, but double precision does provide greater accuracy.  Often the added accuracy wont matter, but I have had research problems where it did, so beware!

===============================================================================
2. Running
-------------------------------------------------------------------------------

There is a sample script in the Run directory that may be of use, but the vast majority of the script is used to configure the job scheduler, and you will be required to tweak that yourself for any system you run on.  The only requirements from the software side is that the code is run from a directory that it can deposit large amounts of data to, and that the location of an appropriate configuration file is handed in as the sole argument to the program.  

The code periodically has three types of outputs, which happens at configurable intervals. The first is a box average of various quantities of interest, such as the peak velocity or the magnetic energy density. See the comments in IO.c for more details on these.  Additionally, it periodically dumps out the full contents of the spatial arrays. Data is laid out as a simple 3D array with the x dimension being contiguous and the z dimension being least contiguous. Finally, the code also has checkpoint outputs where the important sections of memory from each processor is dumped as a list of files. The code keeps around the two most recent dumps, so that even if the code terminates during the writing of a dump, thus corrupting it, a sane restart condition still exists. It should be noted that restarting from a checkpoint is only possible if the same number and layout of processors is used between runs.

The behavior of the code during runtime is determined by a configuration file which must be supplied as the first and only command line argument when the code is launched. An example file is in src/config.cfg. Pairs of [Descriptor] delineate groups of parameters that can be specified, very similar to how Fortran namelists work. Each parameter is specified as a name=value pair. 

===============================================================================
3. Numerical Details
-------------------------------------------------------------------------------

Here is a brief overview of some of the numerical algorithms utilized in this code.  See comments within the code for additional details.

-------------------------------------------------------------------------------
3a. Pseudo-Spectral Algorithms
-------------------------------------------------------------------------------

The code employs a pure pseudo-spectral approach, and solves the equations in a triply periodic cartesian box. The basics of the pseudo-spectral method are as follows. For spectral techniques, each field of interest is stored as a set of Fourier modes, rather than a set of values at points in space. The major advantages of this approach are twofold. First, you gain exponential convergence with the number of wave modes in the problem. Second, derivatives become trivial operations where you multiply the amplitude by the wave-number. Where spectral methods are limited is in the calculation of nonlinear terms, which involve an n^2 number of wave-wave interactions. Pseudo-spectral techniques improve on this by first employing an inverse fft to recover the spatial fields, calculate nonlinear terms here through simply multiplication, and then fft the system back into wave modes. One must be careful with aliasing when using a discrete Fourier transform, but if done correctly the overall effect is to drop the scaling of the nonlinear terms from n^2 to n log n.

-------------------------------------------------------------------------------
3b. Solenoidal Condition
-------------------------------------------------------------------------------

The equations solved by this code have two very important constraints that *must* be observed; namely the divergence of the magnetic field and the velocity field must be exactly zero. To achieve this, each of these vector fields are decomposed into their poloidal and toroidal scalar components, where any divergence free vector field A can be written as A=∇x(T) + ∇x∇x(P) + A(z). These scalar fields are then evolved through time rather than their associated vector fields, and the details of decomposing a vector into these scalar fields can be found in the comments for Numerics.c

-------------------------------------------------------------------------------
3c. Misc Features
-------------------------------------------------------------------------------

-- Time updates with explicit 3rd level Adams-Bashforth method
-- Time step changes dynamically each iteration
-- Each term in the equations can be easily enabled/disabled at runtime
-- Code has a built-in logging system with customizable levels of output (Trace/Debug/Info/Warn/Error)
-- Robust checkpointing system allowing the code to be restarted from previous execution, even if it did not exit cleanly
