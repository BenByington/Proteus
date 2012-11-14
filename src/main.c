/* 
 * File:   main.cpp
 * Author: Ben
 *
 * Created on February 12, 2010, 1:43 PM
 */

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


/*
 * 
 */
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    if(argc != 2)
    {
        fprintf(stderr, "Exactly one command line parameter must be handed in, the location of the configuration file!!\n");
        int i;
        for(i = 0; i < argc; i++)
            fprintf(stderr, "arg %d: %s\n",i, argv[i]);
        return -1;
    }
    srand(time(0));

    setupEnvironment(argv[1]);
    info("Code Initialization Complete\n",0);

    testIO();
    if(compute_node)
        testPT();

    initState();
    initIO();
    if(compute_node)
    {
        initPhysics();
    }

    if(compute_node)
    {
    readSpatial(u->vec->x, "x");
    readSpatial(u->vec->y, "y");
    readSpatial(u->vec->z, "z");

    fftForward(u->vec->x);
    fftForward(u->vec->y);
    fftForward(u->vec->z);

    decomposeSolenoidal(u->sol, u->vec, 0);
    recomposeSolenoidal(u->sol, u->vec);

    fftBackward(u->vec->x);
    fftBackward(u->vec->y);
    fftBackward(u->vec->z);
    
    writeSpatial(u->vec->x, "x2");
    writeSpatial(u->vec->x, "y2");
    writeSpatial(u->vec->x, "z2");
    }
    else
    {
        readSpatial(0, "x");
        readSpatial(0, "y");
        readSpatial(0, "z");
        writeSpatial(0, "x2");
        writeSpatial(0, "y2");
        writeSpatial(0, "z2");
    }

   /*
    while((iteration < maxSteps) && (elapsedTime < maxTime))

    {
        iteration++;
        debug("Working on step %d\n", iteration);
        if(compute_node)
            iterate();

        MPI_Bcast(&elapsedTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        performOutput();
    }*/

    info("Run Complete: Cleaning and Exiting now\n",0);
    if(compute_node)
    {
        finalizePhysics();
        finalizeState();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}

