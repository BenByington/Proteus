/* 
 * File:   Environment.h
 * Author: Ben
 *
 * Created on March 10, 2010, 1:19 PM
 */

#ifndef _ENVIRONMENT_H
#define	_ENVIRONMENT_H

#include <mpi.h>

void setupEnvironment();

extern double PI;

typedef struct
{
    int min;
    int max;
    int width;
}indexes;

//Comms
extern MPI_Comm hcomm;
extern MPI_Comm vcomm;
extern MPI_Comm ccomm;
extern MPI_Comm iocomm;
extern MPI_Comm fcomm;

//layout info
extern int hdiv;
extern int vdiv;
extern int compute_node;

//proc identifications
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

//problem dimensions
extern int nx;
extern int ny;
extern int nz;
extern double xmx;
extern double ymx;
extern double zmx;
extern double dx;
extern double dy;
extern double dz;
extern int nkx;
extern int nky;
extern int nkz;
extern int ndkx;
extern int ndky;
extern int ndkz;

//work distributions
extern indexes * my_x;
extern indexes * all_x;
extern indexes * max_x;
extern indexes * my_z;
extern indexes * all_z;
extern indexes * max_z;
extern indexes * my_kx;
extern indexes * all_kx;
extern indexes * max_kx;
extern indexes dealias_kx;
extern indexes * my_ky;
extern indexes * all_ky;
extern indexes dealias_ky;
extern indexes * max_ky;
extern indexes dealias_kz;
extern int spatialCount;
extern int spectralCount;

//io distributions
extern indexes * io_layers;
extern int my_io_layer;
extern int nz_layers;

//initial conditions
#define SCRATCH 0
#define SPATIAL 1
#define CHECKPOINT 2
extern int startType;
extern char * startDir;

//IO configurations
extern int n_io_nodes;
extern int statusRate;
extern int spatialRate;
extern int scalarRate;
extern int scalarPerF;
extern int checkRate;
extern int checkDir;

//physics terms
extern int momEquation;
extern int magEquation;
extern int tEquation;
extern int momAdvection;
extern int viscosity;
extern int buoyancy;
extern int momStaticForcing;
extern int momTimeForcing;
extern int lorentz;
extern char * forceFile;
extern int tDiff;
extern int tempAdvection;
extern int magDiff;
extern int magAdvect;
extern int kinematic;
extern int magTimeForcing;

//phsyics parameters
extern double Pr;
extern double Ra;
extern double Pm;
extern double alpha;

//integration parameters
extern int maxSteps;
extern double maxTime;
extern int iteration;
extern double safetyFactor;
extern double dt;
extern double dt1;
extern double dt2;
extern double elapsedTime;

#endif	/* _ENVIRONMENT_H */

