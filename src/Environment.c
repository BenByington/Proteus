#include "Environment.h"
#include "LaborDivision.h"
#include "Communication.h"
#include "Log.h"
#include "Properties.h"

#include <stdio.h>
#include <stdlib.h>

void setupEnvironment(char * loc)
{
    PI = 4.0 * atan2(1.0, 1.0);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);

    initLogging();
    info("Comencing Environment Initialization\n",0);

    loadPrefs(loc);

    if(n_io_nodes != gsize - hdiv*vdiv)
    {
        error("ERROR!  Wrong number of processors requested!  n_io_nodes = %d.  The code will now crash gracelessly\n", n_io_nodes);
        abort();
    }

    lab_initGeometry();
    lab_initGroups();
    lab_initDistributions();

    //TODO: Don't forget to turn measuring back on!
    if(compute_node)
    {
        com_init(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

double PI;

MPI_Comm hcomm;
MPI_Comm vcomm;
MPI_Comm ccomm;
MPI_Comm iocomm;
MPI_Comm fcomm;

int hdiv;
int vdiv;
int compute_node;

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
double xmx;
double ymx;
double zmx;
double dx;
double dy;
double dz;
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

int startType;
char * startDir = 0;


int n_io_nodes;
int statusRate = 100;
int spatialRate = 1000;
int scalarRate = 1000;
int scalarPerF = 1;

int forcing = 0;
int momEquation = 0;
int magEquation = 0;
int tEquation = 0;
int momAdvection = 0;
int viscosity = 0;
int buoyancy = 0;
int lorentz = 0;
char * forceFile = 0;
int tDiff = 0;
int tempAdvection = 0;
int magDiff = 0;
int magAdvect = 0;

double Re = 0;
double Ra = 0;
double Prm = 0;
double alpha = 1.0;

int nSteps = 0;
int iteration = 0;
double safetyFactor = 0;
double dt = 0;
double dt1 = 0;
double dt2 = 0;
double elapsedTime = 0;
