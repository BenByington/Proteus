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

#include "Environment.h"
#include "Log.h"
#include <mpi.h>
#include <sys/stat.h>


FILE * procFile;
FILE * progress;

/*
 * Each process gets its own log file in the Logs directory.
 */
void initLogging()
{
    if(grank == 0)
        mkdir("Logs",  S_IRWXU);

    MPI_Barrier(MPI_COMM_WORLD);

    char name[100];
    sprintf(name, "Logs/Proc%d", grank);

    procFile = fopen(name, "w");

    trace("System wide tracing enabled\n");
    debug("System wide debug enabled\n");
    info("System wide info enabled\n");
    warn("System wide warning enabled\n");
    error("System wide error enabled\n");
    
    return;
}

void endLogging()
{
    info("Shutting down log file stream\n");
    fclose(procFile);
    
    return;
}

