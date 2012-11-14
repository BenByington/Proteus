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

int benchmark();
int execute(char * propFile);


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
    
    
    return execute(argv[1]);
}

int execute(char * propLoc)
{
    srand(time(0));

    loadPrefs(propLoc);

    setupEnvironment();
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
    MPI_Finalize();
}

