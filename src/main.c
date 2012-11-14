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

    return status;
}

int execute(char * propLoc)
{
    loadPrefs(propLoc);


    info("Code Initialization Complete\n",0);
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
        debug("Working on step %d\n", iteration);
        if(compute_node)
            iterate();

        MPI_Bcast(&elapsedTime, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);
        performOutput();
    }

    info("Run Complete: Cleaning and Exiting now\n",0);
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

    int contDiv = 1;
    int contSize = 1;
    
    loadPrefs(propLoc);

    hdiv = 3;
    vdiv = 3;
    n_io_nodes = gsize - hdiv * vdiv;
    while(n_io_nodes > 0)
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
                gettimeofday(&stop, NULL);
            }
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
        n_io_nodes = gsize - hdiv * vdiv;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}