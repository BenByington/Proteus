#include "Environment.h"
#include "Log.h"
#include <mpi.h>
#include <sys/stat.h>


FILE * procFile;
FILE * progress;

void initLogging()
{
    if(grank == 0)
        mkdir("Logs",  S_IRWXU);

    MPI_Barrier(MPI_COMM_WORLD);

    char name[100];
    sprintf(name, "Logs/Proc%d", grank);

    procFile = fopen(name, "w");

    trace("System wide tracing enabled\n",0);
    debug("System wide debug enabled\n",0);
    info("System wide info enabled\n",0);
    warn("System wide warning enabled\n",0);
    error("System wide error enabled\n",0);
}

