/*
 * Copywrite 2013 Benjamin Byington
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

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "Field.h"
#include "Communication.h"
#include "Environment.h"
#include "IO.h"
#include "Logs/Log.h"
#include "Numerics.h"
#include "State.h"
#include "Physics.h"
#include "Properties.h"
#include "LaborDivision.h"

int benchmark(char * propFile);
int execute(char * propFile);

#include "cartesian/periodic/FieldPeriodic.h"
#include "cartesian/periodic/VectorPeriodic.h"
#include "cartesian/periodic/SolenoidPeriodic.h"
#include "cartesian/periodic/TensorPeriodic.h"

/*
 * This is the main driving routine.  Doesn't do much save initialize MPI and
 * a few other things, and then call the main execution loop.
 */
int main(int argc, char** argv)
{
    FieldPeriodic temp1();
    VectorPeriodic temp2();
    SolenoidPeriodic temp3();
    TensorPeriodic temp4();
    int status;

    //Start up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);
    
    //Ensure correct calling signature.  Note: Benchmarking is not currently
    //operative.
    if(argc < 2 || argc > 3)
    {
        fprintf(stderr, "Usage: Proteus <propFile> [benchmark]\n");
        int i;
        for(i = 0; i < argc; i++)
            fprintf(stderr, "arg %d: %s\n",i, argv[i]);
        return -1;
    }

    srand(time(0));

    initLogging();
    
    if(argc == 2)
        status = execute(argv[1]);
    else
        status = benchmark(argv[1]);

    //shut it all down.
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    endLogging();
    
    return status;
}

/*
 * This is the main execution loop.  It does not do much besides initializing
 * the physical problem and then entering a tight iteration evolving the system
 * until either the maximum number of iterations has executed or the maximum 
 * amount of simulation time has passed.
 */
int execute(char * propLoc)
{
    loadPrefs(propLoc);

    info("Code Initialization Complete\n");
    setupEnvironment();
    
    //Sanity checks.  These should either be removed or moved to a optional
    //unit test framework.
    testIO();
    if(compute_node)
        testPT();

    initState();
    initIO();
    if(compute_node)
    {
        initPhysics();
    }

    while((iteration < maxSteps) && (elapsedTime < maxTime))
    {
        iteration++;
        info("Working on step %d\n", iteration);
        if(compute_node)
            iterate();

        MPI_Bcast(&elapsedTime, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);
        performOutput();
        
        //This is an experimental section where the domain moves during
        //computation to keep an item of interest centered.  Not fully
        //operational...
//        if(recentering && recenterTerminate != 0)
//        {
//            if((*recenterTerminate)() == 1) break;
//        }
    }

    info("Run Complete: Cleaning and Exiting now\n");
    if(compute_node)
    {
        finalizePhysics();
        finalizeState();
        com_finalize();
    }
    finalizeIO();
    lab_finalize();

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}

/*
 * Not currently functional!
 */
int benchmark(char * propLoc)
{
    double dstart;
    double dstop;
    struct timeval start;
    struct timeval stop;
    
    loadPrefs(propLoc);

    //make sure we don't start running to the file system for data
    startFlag = SCRATCH;
    momStaticForcing = 0;

    hdiv = 3;
    vdiv = 3;
    n_io_nodes = vdiv;
    while(hdiv * vdiv + n_io_nodes <= gsize)
    {
        setupEnvironment();

        initState();
        initIO();
        if(compute_node)
        {
            initPhysics();
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if(compute_node)
        {
            gettimeofday(&start,NULL);
            int i;
            for(i = 0; i < 100; i++)
            {
                fftForward(B->vec->x);
                fftBackward(B->vec->x);
            }
            gettimeofday(&stop, NULL);
            dstart = start.tv_sec+(start.tv_usec/1000000.0);
            dstop = stop.tv_sec + (stop.tv_usec/1000000.0);

            if(grank == 0)
            {
                fprintf(stderr, "hdiv: %d, vdiv: %d, time: %f\n", hdiv, vdiv, dstop - dstart);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /*while((iteration < maxSteps) && (elapsedTime < maxTime))
        {
            iteration++;
            debug("Working on step %d\n", iteration);
            if(compute_node)
                iterate();

            MPI_Bcast(&elapsedTime, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);
            performOutput();
        }*/

        if(compute_node)
        {
            finalizePhysics();
            finalizeState();
            com_finalize();
        }
        lab_finalize();

        hdiv++;
        vdiv++;
        n_io_nodes = vdiv;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}

