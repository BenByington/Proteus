/*
 * Copyright 2013 Benjamin Byington
 *
 * This file is part of the Proteus software package
 * 
 * Proteus is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Proteus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with Proteus.  If not, see <http://www.gnu.org/licenses/>
 */

#include "State.h"
#include "Logs/Log.h"
#include "Environment.h"
#include "IO.h"
#include "Numerics.h"
#include "Communication.h"
#include "Field.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

//These are "private" and never called outside this file.
void startScratch();
void startSpatial();

/*
 * Here we allocate memory for our state variables, and initialize them
 * appropriately.
 */
void initState()
{
    info("Setting up the initial conditions\n");
    if(compute_node)
    {
        /*
         * Note:  The forcing fields never get initialized to 0.  This is safe
         *        because the integration scheme ramps up in accuracy, and wont
         *        try to access any of the "dirty" arrays before they have been
         *        filled with valid data.
         */
        B = newComponentVar();
        u = newComponentVar();
        T = (p_field)malloc(sizeof(field));
        allocateSpatial(T);
        allocateSpectral(T);
        allocateForce(T);
        
        //This is for a hyper diffusion applied to the boundaries to try and
        //zero them out without breaking divergence constraints.  Currently does
        //not work as desired, and is disabled by default.
        if(sanitize)
        {
            hyper = (p_field)malloc(sizeof(field));
            hyperWork = (p_field)malloc(sizeof(field));
            allocateSpectral(hyperWork);
            allocateSpatial(hyperWork);
            allocateSpatial(hyper);
            
            //Use tanh functions to give us function that is 1 close to the
            //boundaries and zero in the interior.  These will be weights 
            //applied to they hyper diffusion.
            int i,j,k;
            int index = 0;
            double temp;
            for(i = my_z->min; i <= my_z->max; i++)
            {
                for(j = my_x->min; j <= my_x->max; j++)
                {
                    if(ny > 1)
                    {
                        for(k = 0; k < ny; k++)
                        {
                            hyper->spatial[index] = .5 * tanh(3 - k) + .5 * tanh(k + 3 - ny) + .5;
                            temp =.5 * tanh(3 - j) + .5 * tanh(j + 3 - nx) + .5;
                            if(temp < hyper->spatial[index])
                                hyper->spatial[index] = temp;
                            temp =.5 * tanh(3 - i) + .5 * tanh(i + 3 - nz) + .5;
                            if(temp < hyper->spatial[index])
                                hyper->spatial[index] = temp;
                            index++;
                        }
                    }
                    else
                    {
                        
                        hyper->spatial[index] = .5 * tanh(3 - j) + .5 * tanh(j + 3 - nx) + .5;
                        temp =.5 * tanh(3 - i) + .5 * tanh(i + 3 - nz) + .5;
                        if(temp < hyper->spatial[index])
                            hyper->spatial[index] = temp;
                        index++;
                    }
                }
            }
        }
    }

    if(startFlag == CHECKPOINT)
    {
        //IO.c knows how to read in data from a checkpoint.
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
            PRECISION ddumb = 0;
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
                //Defined below
                startScratch();
        }
        else if(startFlag == SPATIAL)
        {
            //Defined below.
            startSpatial();
        }
        else
        {
            error("Unsupported start type requested: %d\n", startFlag);
        }
    }
    

    //Time independent forcing for the momentum equation.
    if(momStaticForcing)
    {
        if(!viscosity)
        {
            warn("Static momentum forcing is currently designed to counteract the diffusion term.  Diffusion is disabled so forcing will be neglected for this run\n");
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
                laplacian(forceField->spectral, forceField->spectral, 0, -Pr);
                eraseSpatial(forceField);
            }
        }
    }
    
    //Time independent forcing for the magnetic induction equation.
    if(magStaticForcing)
    {
        if(!magDiff)
        {
            warn("Static magnetic forcing is currently designed to counteract the diffusion term.  Magnetic diffusion is disabled so forcing will be neglected for this run\n");
            magStaticForcing = 0;
        }
        else
        {
            info("Loading forcing file: %s\n", magForceFile);

            if(compute_node)
            {
                magForceField = (p_field)malloc(sizeof(field));
                allocateSpatial(magForceField);
                allocateSpectral(magForceField);
            }

            readSpatial(magForceField, magForceFile);

            if(compute_node)
            {
                fftForward(magForceField);
                laplacian(magForceField->spectral, magForceField->spectral, 0, -Pr / Pm);
                eraseSpatial(magForceField);
            }
        }
    }
}

/*
 * Undo all memory allocations from the above routine.
 */
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
    if(magStaticForcing)
    {
        eraseSpectral(magForceField);
        free(magForceField);
        magForceField = 0;
    }
    
    if(sanitize)
    {
        eraseSpatial(hyper);
        eraseSpatial(hyperWork);
        eraseSpectral(hyperWork);
        free(hyper);
        free(hyperWork);
    }
}

/*
 * If we are starting a simulation from scratch, then all of our state variables
 * need to be initialized to zero at all points in the domain.
 */
void startScratch()
{
    info("Code is starting from scratch\n");

    memset(T->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(T->spectral, 0, spectralCount * sizeof(complex<PRECISION>));

    memset(B->vec->x->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(B->vec->y->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(B->vec->z->spatial, 0, spatialCount * sizeof(PRECISION));

    memset(B->vec->x->spectral, 0, spectralCount * sizeof(complex<PRECISION>));
    memset(B->vec->y->spectral, 0, spectralCount * sizeof(complex<PRECISION>));
    memset(B->vec->z->spectral, 0, spectralCount * sizeof(complex<PRECISION>));

    memset(B->sol->poloidal->spectral, 0, spectralCount * sizeof(complex<PRECISION>));
    memset(B->sol->toroidal->spectral, 0, spectralCount * sizeof(complex<PRECISION>));

    memset(B->sol->mean_x, 0, ndkz * sizeof(complex<PRECISION>));
    memset(B->sol->mean_y, 0, ndkz * sizeof(complex<PRECISION>));

    B->sol->mean_z = 0;

    memset(u->vec->x->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(u->vec->y->spatial, 0, spatialCount * sizeof(PRECISION));
    memset(u->vec->z->spatial, 0, spatialCount * sizeof(PRECISION));

    memset(u->vec->x->spectral, 0, spectralCount * sizeof(complex<PRECISION>));
    memset(u->vec->y->spectral, 0, spectralCount * sizeof(complex<PRECISION>));
    memset(u->vec->z->spectral, 0, spectralCount * sizeof(complex<PRECISION>));

    memset(u->sol->poloidal->spectral, 0, spectralCount * sizeof(complex<PRECISION>));
    memset(u->sol->toroidal->spectral, 0, spectralCount * sizeof(complex<PRECISION>));

    memset(u->sol->mean_x, 0, ndkz * sizeof(complex<PRECISION>));
    memset(u->sol->mean_y, 0, ndkz * sizeof(complex<PRECISION>));

    u->sol->mean_z = 0;
}

/*
 * User has specified a folder on disk containing dumps of the state variables,
 * and we will read them in and use them as initial conditions.
 */
void startSpatial()
{
    char name[100];

    //If we are continuing from another simulation, but not using checkpoints 
    //for some reason, then we need to recover the simulation time that this
    //data dump was made at.
    if(grank == 0)
    {
        sprintf(name,"%s/info",startDir);
        FILE * info;
        info = fopen(name, "r");
        if(info)
        {
            fscanf(info, infostri, &elapsedTime);
            fclose(info);
        }
        else
        {
            elapsedTime = 0;
        }
    }

    //share out the simulation time so everyone knows.
    MPI_Bcast(&elapsedTime, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);

    //Read in state variables for any active equations.  Don't bother for
    //variable that won't be used.  They are read in the spatial coordinates
    //they were stored in, converted to spectral coordinates, and if necessary
    //converted further into toroidal/poloidal decomposition.
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

p_field hyper;
p_field hyperWork;

PRECISION maxVel[3];
p_field forceField = 0;
p_field magForceField = 0;
