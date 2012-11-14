#include "State.h"
#include "Log.h"
#include "Environment.h"
#include "IO.h"
#include "Numerics.h"
#include "Communication.h"
#include "Field.h"

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

    if(startFlag == CHECKPOINT)
    {
        readCheckpoint();
    }
    else
    {
        //do some quick setup so future starts can be a checkpoint
        if(grank == 0)
        {
            //dummy data so the file is the right size, but we can tell
            //nothing has been computed yet and both Checkpoint directories
            //are empty (or have invalid data)
            FILE * out = fopen("Checkpoint0/state", "w");
            int idumb = 0;
            int ddumb = 0;
            fwrite(&ddumb, sizeof(PRECISION), 1, out);
            fwrite(&ddumb, sizeof(PRECISION), 1, out);
            fwrite(&ddumb, sizeof(PRECISION), 1, out);
            fwrite(&idumb, sizeof(int), 1, out);
            fclose(out);

            out = fopen("Checkpoint1/state", "w");
            fwrite(&ddumb, sizeof(PRECISION), 1, out);
            fwrite(&ddumb, sizeof(PRECISION), 1, out);
            fwrite(&ddumb, sizeof(PRECISION), 1, out);
            fwrite(&idumb, sizeof(int), 1, out);
            fclose(out);
        }


        if(startFlag == SCRATCH)
        {
            //no need for IO nodes in here
            if(compute_node)
                startScratch();
        }
        else if(startFlag == SPATIAL)
        {
            startSpatial();
        }
        else
        {
            error("Unsupported start type requested: %d\n", startFlag);
        }
    }
    

    if(momStaticForcing)
    {
        if(!viscosity)
        {
            warn("Static momentum forcing is currently designed to counteract the diffusion term.  Diffusion is disabled so forcing will be neglected for this run\n",0);
            momStaticForcing = 0;
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
                laplacian(forceField->spectral, forceField->spectral, 0, -1.0 / Pr);
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

    if(momStaticForcing)
    {
        eraseSpectral(forceField);
        free(forceField);
        forceField = 0;
    }
}

void startScratch()
{
    info("Code is starting from scratch\n",0);

    memset(T->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(T->spectral, 0, spectralCount * sizeof(complex PRECISION));

    memset(B->vec->x->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(B->vec->y->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(B->vec->z->spatial, 0, spatialCount * sizeof(PRECISION));

    memset(B->vec->x->spectral, 0, spectralCount * sizeof(complex PRECISION));
    memset(B->vec->y->spectral, 0, spectralCount * sizeof(complex PRECISION));
    memset(B->vec->z->spectral, 0, spectralCount * sizeof(complex PRECISION));

    memset(B->sol->poloidal->spectral, 0, spectralCount * sizeof(complex PRECISION));
    memset(B->sol->toroidal->spectral, 0, spectralCount * sizeof(complex PRECISION));

    memset(B->sol->mean_x, 0, ndkz * sizeof(complex PRECISION));
    memset(B->sol->mean_y, 0, ndkz * sizeof(complex PRECISION));

    B->sol->mean_z = 0;

    memset(u->vec->x->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(u->vec->y->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(u->vec->z->spatial, 0, spatialCount * sizeof(PRECISION));

    memset(u->vec->x->spectral, 0, spectralCount * sizeof(complex PRECISION));
    memset(u->vec->y->spectral, 0, spectralCount * sizeof(complex PRECISION));
    memset(u->vec->z->spectral, 0, spectralCount * sizeof(complex PRECISION));

    memset(u->sol->poloidal->spectral, 0, spectralCount * sizeof(complex PRECISION));
    memset(u->sol->toroidal->spectral, 0, spectralCount * sizeof(complex PRECISION));

    memset(u->sol->mean_x, 0, ndkz * sizeof(complex PRECISION));
    memset(u->sol->mean_y, 0, ndkz * sizeof(complex PRECISION));

    u->sol->mean_z = 0;
}

void startSpatial()
{
    char name[100];

    if(compute_node)
    {
        if(magEquation)
        {
            readSpatial(B->vec->x, 0);
            fftForward(B->vec->x);

            readSpatial(B->vec->y, 0);
            fftForward(B->vec->y);

            readSpatial(B->vec->z, 0);
            fftForward(B->vec->z);
            
            decomposeSolenoidal(B->sol, B->vec,0);
        }

        if(momEquation)
        {
            readSpatial(u->vec->x, 0);
            fftForward(u->vec->x);

            readSpatial(u->vec->y, 0);
            fftForward(u->vec->y);

            readSpatial(u->vec->z, 0);
            fftForward(u->vec->z);
            
            decomposeSolenoidal(u->sol, u->vec,0);
        }

        if(tEquation)
        {
            readSpatial(T, 0);
            fftForward(T);
        }
        
    }
    else if(io_node)
    {
        if(magEquation)
        {
            sprintf(name,"%s/Bx",startDir);
            readSpatial(0, name);

            sprintf(name,"%s/By",startDir);
            readSpatial(0, name);

            sprintf(name,"%s/Bz",startDir);
            readSpatial(0, name);
        }

        if(momEquation)
        {
            sprintf(name,"%s/u",startDir);
            readSpatial(0, name);

            sprintf(name,"%s/v",startDir);
            readSpatial(0, name);

            sprintf(name,"%s/w",startDir);
            readSpatial(0, name);
        }

        if(tEquation)
        {
            sprintf(name,"%s/T",startDir);
            readSpatial(0, name);
        }
    }

}

p_componentVar B;
p_componentVar u;
p_field T;

PRECISION maxVel[3];
p_field forceField = 0;

