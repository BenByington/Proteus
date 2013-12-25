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

/* 
 * File:   main.cpp
 * Author: Ben
 *
 * Created on February 12, 2010, 1:43 PM
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
#include "Log.h"
#include "Numerics.h"
#include "State.h"
#include "Physics.h"
#include "Properties.h"
#include "LaborDivision.h"

int benchmark(char * propFile);
int execute(char * propFile);


/*
 * 
 */
int main(int argc, char** argv)
{
    int status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);
    
    if(argc < 2 || argc > 3)
    {
        fprintf(stderr, "Usage: imhd <propFile> [benchmark]\n");
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

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    endLogging();
    
    return status;
}

int execute(char * propLoc)
{
    loadPrefs(propLoc);

    info("Code Initialization Complete\n");
    setupEnvironment();
    
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
    lab_finalize();

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}

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

