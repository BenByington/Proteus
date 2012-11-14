#include "State.h"
#include "Log.h"
#include "Environment.h"
#include "IO.h"
#include "Numerics.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

void startScratch();
void startSpatial();

void initState()
{
    info("Setting up the initial conditions\n",0);
    if(compute_node)
    {
        B = newComponentVar();
        u = newComponentVar();
        T = (p_field)malloc(sizeof(field));
        allocateSpatial(T);
        allocateSpectral(T);
        allocateForce(T);
    }

    if(startType == SCRATCH)
    {
        //no need for IO nodes in here
        if(compute_node)
            startScratch();
    }
    else if(startType == SPATIAL)
    {
        startSpatial();
    }
    else
    {
        error("Unsupported start type requested: %d\n", startType);
    }
    

    if(forcing)
    {
        if(!viscosity)
        {
            warn("Forcing is currently designed to counteract the diffusion term.  Diffusion is disabled so forcing will be neglected for this run\n",0);
            forcing = 0;
        }
        else
        {
            info("Loading forcing file: %s\n", forceFile);

            if(compute_node)
            {
                forceField = (p_field)malloc(sizeof(field));
                allocateSpatial(forceField);
                allocateSpectral(forceField);
            }

            readSpatial(forceField, forceFile);

            if(compute_node)
            {
                fftForward(forceField);
                laplacian(forceField->spectral, forceField->spectral, 0, -1.0 / Re);
                eraseSpatial(forceField);
            }
        }
    }
}

void finalizeState()
{
    deleteComponentVar(&B);
    deleteComponentVar(&u);
    eraseSpatial(T);
    eraseSpectral(T);
    eraseForce(T);
    free(T);

    if(forcing)
    {
        eraseSpectral(forceField);
        free(forceField);
        forceField = 0;
    }
}

void startScratch()
{
    info("Code is starting from scratch\n",0);

    memset(T->spatial, 0, spatialCount * sizeof(double));
    memset(T->spectral, 0, spectralCount * sizeof(complex double));

    memset(B->vec->x->spatial, 0, spatialCount * sizeof(double));
    memset(B->vec->y->spatial, 0, spatialCount * sizeof(double));
    memset(B->vec->z->spatial, 0, spatialCount * sizeof(double));

    memset(B->vec->x->spectral, 0, spectralCount * sizeof(complex double));
    memset(B->vec->y->spectral, 0, spectralCount * sizeof(complex double));
    memset(B->vec->z->spectral, 0, spectralCount * sizeof(complex double));

    memset(B->sol->poloidal->spectral, 0, spectralCount * sizeof(complex double));
    memset(B->sol->toroidal->spectral, 0, spectralCount * sizeof(complex double));

    memset(B->sol->mean_x, 0, ndkz * sizeof(complex double));
    memset(B->sol->mean_y, 0, ndkz * sizeof(complex double));

    B->sol->mean_z = 0;

    memset(u->vec->x->spatial, 0, spatialCount * sizeof(double));
    memset(u->vec->y->spatial, 0, spatialCount * sizeof(double));
    memset(u->vec->z->spatial, 0, spatialCount * sizeof(double));

    memset(u->vec->x->spectral, 0, spectralCount * sizeof(complex double));
    memset(u->vec->y->spectral, 0, spectralCount * sizeof(complex double));
    memset(u->vec->z->spectral, 0, spectralCount * sizeof(complex double));

    memset(u->sol->poloidal->spectral, 0, spectralCount * sizeof(complex double));
    memset(u->sol->toroidal->spectral, 0, spectralCount * sizeof(complex double));

    memset(u->sol->mean_x, 0, ndkz * sizeof(complex double));
    memset(u->sol->mean_y, 0, ndkz * sizeof(complex double));

    u->sol->mean_z = 0;
}

void startSpatial()
{
    char name[100];

    if(compute_node)
    {
        readSpatial(B->vec->x, 0);
        fftForward(B->vec->x);

        readSpatial(B->vec->y, 0);
        fftForward(B->vec->y);

        readSpatial(B->vec->z, 0);
        fftForward(B->vec->z);

        readSpatial(u->vec->x, 0);
        fftForward(u->vec->x);

        readSpatial(u->vec->y, 0);
        fftForward(u->vec->y);

        readSpatial(u->vec->z, 0);
        fftForward(u->vec->z);

        readSpatial(T, 0);
        fftForward(T);

        decomposeSolenoidal(B->sol, B->vec,0);
        decomposeSolenoidal(u->sol, u->vec,0);
    }
    else
    {
        sprintf(name,"%s/Bx",startDir);
        readSpatial(0, name);

        sprintf(name,"%s/By",startDir);
        readSpatial(0, name);

        sprintf(name,"%s/Bz",startDir);
        readSpatial(0, name);

        sprintf(name,"%s/u",startDir);
        readSpatial(0, name);

        sprintf(name,"%s/v",startDir);
        readSpatial(0, name);

        sprintf(name,"%s/w",startDir);
        readSpatial(0, name);

        sprintf(name,"%s/T",startDir);
        readSpatial(0, name);
    }

}

p_componentVar B;
p_componentVar u;
p_field T;

double maxVel[3];
p_field forceField = 0;

