/*
 * Copywrite 2013 Benjamin Byington
 *
 * This file is part of the IMHD software package
 * 
 * IMHD is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public Liscence as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * IMHD is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with IMHD.  If not, see <http://www.gnu.org/licenses/>
 */

/*********
 * These are all of the global variables used in the program.  Some of these
 * variables are initialized when the configuration file is parsed, and the rest
 * are set during an initialization routine also declared here.
 *********/

#ifndef _ENVIRONMENT_H
#define	_ENVIRONMENT_H

#include "Precision.h"

#include <mpi.h>

/*
 * This routine begins initialization of many of the variables declared here.
 * This must be called exactly once, after the configuration file describing
 * the problem has been parsed and before any real computations begin.
 */
void setupEnvironment();

extern PRECISION PI;

//This helps map between data stored locally and the logical arrays distributed
//across different processes.
typedef struct
{
    int min;
    int max;
    int width;
}indexes;

//describes a displacement vector.  Only used in an experimental code feature
//which may eventually get removed...
typedef struct
{
    PRECISION dx;
    PRECISION dy;
    PRECISION dz;
}displacement;

//We have several different communication groups outside the standard global
//group.  Processor layout is conceptually thought of as a 2D grid of compute 
//nodes along with a 1D group of IO nodes that each own integer layers of 
//compute nodes.
extern MPI_Comm hcomm;  //communication within a given row of compute nodes
extern MPI_Comm vcomm;  //communication within a given column of compute nodes
extern MPI_Comm ccomm;  //communication within the set of all compute nodes
extern MPI_Comm iocomm; //communication with an io node and its compute layer(s)
extern MPI_Comm fcomm;  //Communication between all io nodes

//layout info
extern int hdiv;         //number of compute nodes in a row
extern int vdiv;         //number of compute nodes in a column
extern int compute_node; //number of compute nodes (hdiv * vdiv)
extern int io_node;      //number of io nodes

//Rank and size identifiers for each communication groups a processor may 
//belong to.  The initial letter matches with the initial letter of the 
//corresponding MPI::Comm variable.
extern int grank;
extern int gsize;
extern int hrank;
extern int hsize;
extern int vrank;
extern int vsize;
extern int csize;
extern int crank;
extern int frank;
extern int fsize;
extern int iorank;
extern int iosize;

//problem resolution in gridpoints
extern int nx;
extern int ny;
extern int nz;

//physical extent of the domain
extern PRECISION xmx;
extern PRECISION ymx;
extern PRECISION zmx;

//physical extent of each grid cell
extern PRECISION dx;
extern PRECISION dy;
extern PRECISION dz;

//number of wave modes before de-aliasing
extern int nkx;
extern int nky;
extern int nkz;

//number of wave modes after de-aliasing (see LaborDivision.h)
extern int ndkx;
extern int ndky;
extern int ndkz;

//Variables describing the domain decomposition.  In spatial corrdinates, x and 
//z are parallel.  In spectral coordinates, kx and ky are parallel.  For each 
//of these, we have an indexes structs with the suffic 'my', a'll' and 'max'.  
//'my' describes where this process fits in the logical global array, 'all'
//has the same information for each processor, contained in an array indexed
//by processor rank, and 'max' can be used to find the largest subdomain kept
//by any processor (in case things don't divide evenly)
extern indexes * my_x;
extern indexes * all_x;
extern indexes * max_x;

extern indexes * my_z;
extern indexes * all_z;
extern indexes * max_z;

extern indexes * my_kx;
extern indexes * all_kx;
extern indexes * max_kx;

extern indexes * my_ky;
extern indexes * all_ky;
extern indexes * max_ky;

//describes which wavenumbers are to be dealised after an FFT
extern indexes dealias_kx;
extern indexes dealias_ky;
extern indexes dealias_kz;

//Describes LOCAL size for spatial and spectral arrays.  
//i.e. spatialCount != nx*ny*nz
extern int spatialCount;
extern int spectralCount;

//Describes how the computational grid is distributed among IO nodes
extern indexes * io_layers;
extern int my_io_layer;
extern int nz_layers;

//initial conditions
#define SCRATCH 0          //starting with all fields initialized to 0
#define SPATIAL 1          //starting from specified spatial output
#define CHECKPOINT 2       //starting from latest checkpoint dump
extern char * startType;   
extern int startFlag;
extern char * startDir;    //Directory containing IC if starting from Spatial

//IO configurations
extern int n_io_nodes;
extern int statusRate;     //iterations between updating status file
extern int spatialRate;    //iteration between spatial dumps
extern int scalarRate;     //iterations between scalar reductions
extern int scalarPerF;     //number of scalar outputs to be placed in one file
extern int checkRate;      //How frequently to save simulation state
extern int checkDir;       //Checkpointing alternates between two directions.

//physics terms
extern int momEquation;
extern int magEquation;
extern int tEquation;
extern int momAdvection;
extern int viscosity;
extern int buoyancy;
extern int magBuoy;
extern int lorentz;
extern int tDiff;
extern int tempAdvection;
extern int tempBackground;
extern int magDiff;
extern int magAdvect;

//forcing terms
extern int momStaticForcing;
extern int magStaticForcing;
extern char * forceFile;
extern char * magForceFile;
extern int kinematic;
extern int momTimeForcing;
extern int magTimeForcing;
extern PRECISION momOmega;
extern PRECISION momEps;
extern PRECISION magK;
extern PRECISION magW;
extern PRECISION magB0;

//domain shifting parameters.  Experimental!!!
extern int recentering;
extern int sanitize;
#define NOCENTERING 0
#define BYMAXCENTER 1
#define NOTERMINATE 0
extern int recenterTerminate;

//phsyics parameters
extern PRECISION Pr;
extern PRECISION Ra;
extern PRECISION Pm;
extern PRECISION alpha;
extern PRECISION magBuoyScale;

//integration parameters
extern int maxSteps;              //end simulation after this many iterations
extern PRECISION maxTime;         //end simulation after this much sim time
extern int iteration;             
extern PRECISION safetyFactor;
extern PRECISION dt;
extern PRECISION dt1;
extern PRECISION dt2;
extern PRECISION elapsedTime;

extern char infostro[];
extern char infostri[];

#endif	/* _ENVIRONMENT_H */

