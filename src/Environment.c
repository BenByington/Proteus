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

#include "Environment.h"
#include "LaborDivision.h"
#include "Communication.h"
#include "Log.h"
#include "Properties.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

/*
 * Here we take care of initialization things.  Namely, we set the variables
 * describing the physical extent of the problem, set the variables controlling
 * communication between processes, set the variables controlling the
 * distribution of work, and finally ensure the existence of the directories
 * this program expects to exist.
 */
void setupEnvironment()
{
    PI = 4.0 * atan2(1.0, 1.0);

    if(n_io_nodes > gsize - hdiv*vdiv)
    {
        error("ERROR!  Too many io nodes requested!  n_io_nodes = %d.  The code will now crash gracelessly\n", n_io_nodes);
        abort();
    }

    iteration = 0;
    elapsedTime = 0;

    //Check out LaborDivision.c for most of the initialization code.
    lab_initGeometry();
    lab_initGroups();
    lab_initDistributions();

    if(grank == 0)
    {
        mkdir("Spatial", S_IRWXU);
        mkdir("Scalars", S_IRWXU);
        mkdir("Checkpoint0", S_IRWXU);
        mkdir("Checkpoint1", S_IRWXU);
    }

    //TODO: FFTW speed measuring is turned off.  Make this a configurable
    //parameter!
    if(compute_node)
    {
        com_init(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

PRECISION PI;

MPI_Comm hcomm;
MPI_Comm vcomm;
MPI_Comm ccomm;
MPI_Comm iocomm;
MPI_Comm fcomm;

int hdiv;
int vdiv;
int compute_node;
int io_node;

int grank = -1;
int gsize = -1;
int hrank = -1;
int hsize = -1;
int vrank = -1;
int vsize = -1;
int csize = -1;
int crank = -1;
int frank = -1;
int fsize = -1;
int iorank = -1;
int iosize = -1;

int nx;
int ny;
int nz;
PRECISION xmx;
PRECISION ymx;
PRECISION zmx;
PRECISION dx;
PRECISION dy;
PRECISION dz;
int nkx;
int nky;
int nkz;
int ndkx;
int ndky;
int ndkz;

indexes * my_x;
indexes * all_x;
indexes * max_x;
indexes * my_z;
indexes * all_z;
indexes * max_z;
indexes * my_kx;
indexes * all_kx;
indexes * max_kx;
indexes dealias_kx;
indexes * my_ky;
indexes * all_ky;
indexes dealias_ky;
indexes * max_ky;
indexes dealias_kz;
int spatialCount;
int spectralCount;

indexes * io_layers;
int my_io_layer;
int nz_layers;

char * startType = 0;
int startFlag = SCRATCH;
char * startDir = 0;


int n_io_nodes;
int statusRate = 100;
int spatialRate = 1000;
int scalarRate = 1000;
int scalarPerF = 1;
int checkRate = 1000;
int checkDir = 0;

int momEquation = 0;
int magEquation = 0;
int tEquation = 0;
int momAdvection = 0;
int viscosity = 0;
int buoyancy = 0;
int magBuoy = 0;
int magBuoyTemp = 1;  //similar to tempBackground below.  Should be on for consistency, but can be disabled
int lorentz = 0;
int tDiff = 0;
int tempAdvection = 0;
int tempBackground = 1; //this controls an if statement nested in an if reliant on the tempAdvection flag.  We want to default them to act together, as this is the norm.
int magDiff = 0;
int magAdvect = 0;

//forcing terms
int momStaticForcing = 0;
int magStaticForcing = 0;
char * forceFile = 0;
char * magForceFile = 0;
int kinematic = 0;
int momTimeForcing = 0;
int magTimeForcing = 0;
PRECISION momOmega = 0;
PRECISION momEps = 0;
PRECISION magK = 0;
PRECISION magW = 0;
PRECISION magB0 = 0;

PRECISION Pr = 0;
PRECISION Ra = 0;
PRECISION Pm = 0;
PRECISION alpha = 1.0;
PRECISION magBuoyScale = 1.0;

//domain shifting parameters
int recentering = NOCENTERING;
int sanitize = 0;
int recenterTerminate = NOTERMINATE;

int maxSteps = 0;
PRECISION maxTime = 0;
int iteration = 0;
PRECISION safetyFactor = 0;
PRECISION dt = 0;
PRECISION dt1 = 0;
PRECISION dt2 = 0;
PRECISION elapsedTime = 0;

char infostro[] = "Time: %g\n";
#ifdef FP
char infostri[] = "Time: %f\n";
#else
char infostri[] = "Time: %lf\n";
#endif

