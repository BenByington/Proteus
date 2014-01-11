/*
 * Copywrite 2013 Benjamin Byington
 *
 * This file is part of the IMHD software package
 * 
 * IMHD is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include "IO.h"
#include "Field.h"
#include "Environment.h"
#include "Logs/Log.h"
#include "State.h"
#include "Numerics.h"
#include "Communication.h"

using namespace std;

FILE * status = 0;

int scalarCount;
PRECISION * scalarData;
PRECISION * piScalarData;
int numScalar;

/*
 * This is a very rudimentary test routine to ensure that we can write data
 * to disk and read it back in without corruption.  We just create three simple
 * data fields with simple spatial dependencies, and write/read them.
 */
void testIO()
{
    int i,j,k;
    field * f = (field *)malloc(3 * sizeof(field));
    field * f2 = (field *)malloc(3 * sizeof(field));

    info("Testing IO routines\n");
    if(grank == 0)
        mkdir("Test", S_IRWXU);

    if(compute_node)
    {
        trace("Creating data for test IO\n");
        allocateSpatial(f);
        allocateSpatial(f+1);
        allocateSpatial(f+2);
        allocateSpatial(f2);
        allocateSpatial(f2+1);
        allocateSpatial(f2+2);

        //f[0](x,y,z) = x
        //f[1](x,y,z) = y
        //f[2](x,y,z) = z
        for(i = 0; i < my_z->width; i++)
        {
            for(j = 0; j < my_x->width; j++)
            {
                for(k = 0; k < ny; k++)
                {
                    int index = k + j*ny + i*my_x->width*ny;
                    f[0].spatial[index] = j + my_x->min;
                    f[1].spatial[index] = k;
                    f[2].spatial[index] = i + my_z->min;
                }
            }
        }

    }

    debug("Writing test IO data to files\n");
    char * name1 = new char[15];
    char * name2 = new char[15];
    char * name3 = new char[15];
    strcpy(name1, "Test/x");
    strcpy(name2, "Test/y");
    strcpy(name3, "Test/z");
    writeSpatial(f, name1);
    writeSpatial(f+1, name2);
    writeSpatial(f+2, name3);

    debug("Reading test IO data from files\n");
    readSpatial(f2, name1);
    readSpatial(f2+1, name2);
    readSpatial(f2+2, name3);

    delete [] name1;
    delete [] name2;
    delete [] name3;
    
    if(compute_node)
    {
        debug("Checking integrity of data\n");
        for(i = 0; i < my_z->width; i++)
        {
            for(j = 0; j < my_x->width; j++)
            {
                for(k = 0; k < ny; k++)
                {
                    int index = k + j*ny + i*my_x->width*ny;
                    if(!(f2[0].spatial[index] == j + my_x->min))
                    {
                        error("IO for x failed %g %d!\n",f2[0].spatial[index], j + my_x->min);
                        abort();
                    }
                    if(!(f2[1].spatial[index] == k))
                    {
                        error("IO for y failed %g %d!\n",f2[1].spatial[index], k);
                        abort();
                    }
                    if(!(f2[2].spatial[index] == i + my_z->min))
                    {
                        error("IO for z failed %g %d!\n",f2[2].spatial[index],i + my_z->min);
                        abort();
                    }
                }
            }
        }

        eraseSpatial(f);
        eraseSpatial(f+1);
        eraseSpatial(f+2);
        eraseSpatial(f2);
        eraseSpatial(f2+1);
        eraseSpatial(f2+2);
    }
    free(f);
    free(f2);
    info("IO Test complete\n");
}

/*
 * There are three stages of execution in this routine.  
 * 
 * 1.   Data in an IO group is gathered together.  An IO group consists of an
 *      integer number of compute node layers and a single IO node, and the IO
 *      node is where the data is gathered.
 * 
 * 2.   The data is transposed to have the desired layout in memory.  Data on
 *      a compute node is stored as [z][x][y].  After gathering to the IO node,
 *      the ordering of compute nodes results in a data layout of 
 *      [l][h][z][y][x] where l iterates over layers and h iterates over rows.
 *      We wish to transpose it to be [z][y][x], both because this is a more
 *      intuitive ordering for data during analysis, and because z is the final
 *      remaining distributed dimension, and a raw combination of the data in
 *      all the IO nodes will now result in a well ordered layout.
 * 
 * 3.   The set of all IO nodes perform a parallel write to disk, resulting in
 *      a single file with an expected ordering.
 * 
 *      Note:  The reason compute nodes store data as [z][x][y] instead of 
 *             [z][y][x] is so that after an FFT operation (and it's required
 *             transposed) spectral modes are stored in [kz][ky][kx].  This 
 *             is a simpler ordering to remember, and in the code we are far
 *             more likely to iterate over individual dimensions in spectral
 *             coordinates than spatial ones. IO happens rarely enough that
 *             the cost of the extra transpose required here is probably
 *             negligible, though this should be verified.
 */
void writeSpatial(field * f, char name[])
{
    int i,j,k,l,m;
    debug("Writing spatial data to file %s\n", name);

    int sndcnt = 0;
    PRECISION * rcvbuff = 0;
    PRECISION * sndbuff = 0;
    
    //the extra +1 just gives a little extra room to do an extra loop below.
    //The extra element means nothing, it just makes the code a hair easier
    //to write
    int displs[iosize+1];
    int rcvcounts[iosize];

    debug("consolidating data to IO nodes\n");
    if(compute_node)
    {
        sndcnt = my_x->width * my_z->width * ny;
        trace("Sending %d PRECISIONs\n", sndcnt);
        MPI_Gatherv(f->spatial, sndcnt, MPI_PRECISION, 0, 0, 0, MPI_PRECISION, 0, iocomm);
        debug("Write Spatial completed\n");
        return;
    }
    else if(io_node)
    {
        rcvbuff = (PRECISION *)malloc(nx * ny * nz_layers * sizeof(PRECISION));
        sndbuff = (PRECISION *)malloc(nx * ny * nz_layers * sizeof(PRECISION));
        trace("Total local data will be %d PRECISIONs\n", nx*ny*nz_layers);

        //We need to calculate the starting index that data from each compute
        //processor will begin at in our array.
        //Note:  Since our own IO node is not contributing any data, both our
        //       IO node and the first compute node get to start at a 
        //       displacement of 0.
        displs[0] = 0;
        displs[1] = 0;
        rcvcounts[0] = 0;

        //staggered loop.  We calculate how much data we receive from one
        //processor at the same time we calculate where the data for the
        //next processor will begin storage.
        int * pidspls = displs + 2;
        int * pircvcounts = rcvcounts+1;
        for(i = io_layers[my_io_layer].min; i <= io_layers[my_io_layer].max; i++)
        {
            for(j = 0; j < hdiv; j++)
            {
                *pircvcounts = all_x[j].width * all_z[i].width * ny;
                *pidspls = *(pidspls-1) + *pircvcounts;
                trace("Proc %d should send %d PRECISIONs at displacement %d\n", hdiv * i + j, *pircvcounts, *pidspls);
                pidspls++;
                pircvcounts++;
            }
        }
        MPI_Gatherv(0, 0, MPI_PRECISION, rcvbuff, rcvcounts, displs, MPI_PRECISION, 0, iocomm);
    }

    debug("transposing data so it is properly contiguous\n");
    //rcvbuff is [l][h][vz][hx][y]
    //we want [lz][y][x]
    int indexr = 0;
    int indexs = 0;
    for(i = 0; i < io_layers[my_io_layer].width; i++)
    {
        for(j = 0; j < hdiv; j++)
        {
            int vz = all_z[i + io_layers[my_io_layer].min].width;
            int vzmin = all_z[i + io_layers[my_io_layer].min].min;
            int vzstart = all_z[io_layers[my_io_layer].min].min;
            for(k = 0; k < vz; k++)
            {
                int hx = all_x[j].width;
                int hxmin = all_x[j].min;
                for(l = 0; l < hx; l++)
                {
                    for(m = 0; m < ny; m++)
                    {
                        indexs = ((k + vzmin - vzstart)*ny + m)*nx + l + hxmin;
                        sndbuff[indexs] = rcvbuff[indexr];
                        indexr++;
                    }
                }
            }
        }
    }



    debug("Performing parallel file write\n");
    //TODO: revisit MPI_MODE_SEQUENTIAL and MPI_INFO_NULL to make sure these are what we want
        MPI_File fh;
    MPI_File_open(fcomm, name, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    debug("MPI File opened successfully\n");
    
    //Calculate displacements for each IO processor into the full file.
    int disp = 0;
    
    //loop over each IO processor before us
    for(i = 0; i < my_io_layer; i++)
    {
        //calculate how many layers those IO processors contribute
        for(j = io_layers[i].min; j <= io_layers[i].max; j++)
        {
            disp += all_z[j].width;
        }
    }
    //Convert layers into actual data size.
    disp *= nx * ny * sizeof(PRECISION);
    
    trace("Our view starts at element %d\n", disp);
    trace("Setting view...\n");
    char * type = new char[20];
    strcpy(type, "native");
    MPI_File_set_view(fh, disp, MPI_PRECISION, MPI_PRECISION, type, MPI_INFO_NULL);
    delete [] type;
    trace("Writing to file...\n");
    MPI_File_write(fh, sndbuff, nx * ny * nz_layers, MPI_PRECISION, MPI_STATUS_IGNORE );
    MPI_File_close(&fh);

    free(sndbuff);
    free(rcvbuff);

    debug("Write Spatial completed\n");

}

/*
 * This function is just the inverse of writeSpatial.  See comments for above
 * function.
 */
void readSpatial(field * f, char name[])
{
    int i,j,k,l,m;
    debug("Reading spatial data from file %s\n", name);

    PRECISION * sndbuff = 0;
    PRECISION * rcvbuff = 0;

    if(io_node)
    {
        sndbuff = (PRECISION *)malloc(nx * ny * nz_layers * sizeof(PRECISION));
        rcvbuff = (PRECISION *)malloc(nx * ny * nz_layers * sizeof(PRECISION));
        trace("Total local data will be %d PRECISIONs\n", nx*ny*nz_layers);

        //do a mpi IO operation
        //TODO: revisit MPI_MODE_SEQUENTIAL and MPI_INFO_NULL to make sure these are what we want
        MPI_File fh;
        debug("Reading file\n");
        MPI_File_open(fcomm, name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

        int disp = 0;
        for(i = 0; i < my_io_layer; i++)
        {
            for(j = io_layers[i].min; j <= io_layers[i].max; j++)
            {
                disp += all_z[j].width;
            }
        }
        disp *= nx * ny * sizeof(PRECISION);
        trace("Our file view starts at displacement %d\n", disp);
        trace("Setting view\n"); 
        char * type = new char[20];
        strcpy(type, "native");
        MPI_File_set_view(fh, disp, MPI_PRECISION, MPI_PRECISION, type, MPI_INFO_NULL);
        delete [] type;
        trace("Reading file\n");
        MPI_File_read(fh, sndbuff, nx * ny * nz_layers, MPI_PRECISION, MPI_STATUS_IGNORE );

        MPI_File_close(&fh);

        
        debug("transposing the data for scatter to compute nodes\n");
        //rcvbuff is [l][h][vz][hx][y]
        //sndbuff [l][vz][y][h][hx]
        int indexr = 0;
        int indexs = 0;
        for(i = 0; i < io_layers[my_io_layer].width; i++)
        {
            for(j = 0; j < hdiv; j++)
            {
                int vz = all_z[i + io_layers[my_io_layer].min].width;
                int vzmin = all_z[i + io_layers[my_io_layer].min].min;
                int vzstart = all_z[io_layers[my_io_layer].min].min;
                for(k = 0; k < vz; k++)
                {
                    int hx = all_x[j].width;
                    int hxmin = all_x[j].min;
                    for(l = 0; l < hx; l++)
                    {
                        for(m = 0; m < ny; m++)
                        {
                            indexs = ((k + vzmin - vzstart)*ny + m)*nx + l + hxmin;

                            rcvbuff[indexr] = sndbuff[indexs];
                            indexr++;
                        }
                    }
                }
            }
        }
    }

    int rcvcnt;
    int displs[iosize+1];
    int sndcounts[iosize];
    debug("scattering data to the compute nodes\n");
    if(compute_node)
    {
        rcvcnt = my_x->width * my_z->width * ny;
        trace("I expect to receive %d PRECISIONs\n", rcvcnt);
        MPI_Scatterv(0, 0, 0, MPI_PRECISION, f->spatial, rcvcnt, MPI_PRECISION, 0, iocomm);
    }
    else if(io_node)
    {
        displs[0] = 0;
        displs[1] = 0;
        sndcounts[0] = 0;

        int * pidspls = displs + 2;
        int * pisndcounts = sndcounts+1;
        for(i = io_layers[my_io_layer].min; i <= io_layers[my_io_layer].max; i++)
        {
            for(j = 0; j < hdiv; j++)
            {
                *pisndcounts = all_x[j].width * all_z[i].width * ny;
                *pidspls = *(pidspls-1) + *pisndcounts;
                trace("Sending %d PRECISIONs from displacement %d to proc %d\n", *pisndcounts, *pidspls, i*j);
                pidspls++;
                pisndcounts++;
            }
        }
        MPI_Scatterv(rcvbuff, sndcounts, displs, MPI_PRECISION, 0, 0, MPI_PRECISION, 0, iocomm);
    }

    if(io_node)
    {
        free(rcvbuff);
        free(sndbuff);
    }
    debug("Reading from file done\n");
}

/*
 * There really are only two things we need to initialize.
 * 
 * 1.   How many scalars there will be, which depends on if magnetic fields
 *      are included in the calculation or not.
 * 2.   Set up the array that will hold the scalar values between outputs to
 *      disk.
 */
void initIO()
{
    if(magEquation)
    {
        numScalar = 24;
    }
    else
    {
        numScalar = 13;
    }

    if(crank == 0)
    {
        if(startFlag != CHECKPOINT)
        {
            status = fopen("status", "w");
            fprintf(status, "Stupid message here and now to make me put a nice and informative one later\n\n");
            fclose(status);
        }
        else
        {
            status = fopen("status", "a");
            fprintf(status, "\nRestarting from most recent Checkpoint\n\n");
            fclose(status);
        }

        scalarCount = 0;

        scalarData = (PRECISION*)malloc(numScalar * scalarPerF * sizeof(PRECISION));
        piScalarData = scalarData;
    }
}

/*
 * Clean up the memory that the init method allocated.
 */
void finalizeIO()
{
    if(crank == 0)
    {
        free(scalarData);
    }
}

/*
 * This is the basic entry point for most IO operations.  It should be called
 * once per iteration of the code, and it will determine which types of IO
 * IO should be performed.  The different types of IO are:
 * 
 * 1.   Spatial dumps
 * 2.   Checkpoint dumps
 * 3.   Scalar reductions
 * 4.   Status file updates
 * 
 * The frequency of each of these is controlled by parameters read in from the
 * configuration file.
 */
void performOutput()
{

    if(crank == 0)
    {
        //update the status file?
        if(iteration % statusRate == 0)
        {
            status = fopen("status", "a");
            fprintf(status, "Iteration %d:\n dt = %g\tElapsed Time = %g\nMax Vel: %g %g %g\n", iteration, dt, elapsedTime, maxVel[0], maxVel[1], maxVel[2]);
            fclose(status);
        }
    }

    /*
     * The scalars are single values that result from a global operation on the
     * domain.  The scalars are:
     * 
     * [0]  --  iteration
     * [1]  --  elapsed time
     * [2]  --  min  u
     * [3]  --  max  u
     * [4]  --  mean u
     * [5]  --  min  v
     * [6]  --  max  v
     * [7]  --  mean v
     * [8]  --  min  w
     * [9]  --  max  w
     * [10] --  mean w
     * [11] --  mean kinetic energy density
     * [12] --  max kinetic energy density
     * 
     * ------------ Only present if magnetic fields are included:
     * [13] --  min  Bx
     * [14] --  max  Bx
     * [15] --  mean Bx
     * [16] --  min  By
     * [17] --  max  By
     * [18] --  mean By
     * [19] --  min  Bz
     * [20] --  max  Bz
     * [21] --  mean Bz
     * [22] --  mean magnetic energy density
     * [23] --  max  magnetic energy density* 
     */
    if(compute_node)
    {
        if(iteration % scalarRate == 0)
        {
            PRECISION * local = (PRECISION*)malloc(sizeof(PRECISION)*numScalar);
            int i;

            PRECISION * datax = u->vec->x->spatial;
            PRECISION * datay = u->vec->y->spatial;
            PRECISION * dataz = u->vec->z->spatial;
            PRECISION temp;

            //I feel like there should be an MPI operation that doesn't require
            //me to manually take a local min/max/etc, but I'm failing to find
            //it.  The general MPI_Reduce does an element by element reduction,
            //which is not what we want.
            
            local[0] = iteration;
            local[1] = elapsedTime;
            local[2] = datax[0];
            local[3] = datax[0];
            local[4] = datax[0];
            local[5] = datay[0];
            local[6] = datay[0];
            local[7] = datay[0];
            local[8] = dataz[0];
            local[9] = dataz[0];
            local[10] = dataz[0];
            temp = pow(datax[0],2) + pow(datay[0],2) + pow(dataz[0],2);
            local[11] = temp;
            local[12] = temp;
            for(i = 1; i < spatialCount; i++)
            {
                local[2] = fmin(local[2], datax[i]);
                local[3] = fmax(local[3], datax[i]);
                local[4] += datax[i];
                local[5] = fmin(local[5], datay[i]);
                local[6] = fmax(local[6], datay[i]);
                local[7] += datay[i];
                local[8] = fmin(local[8], dataz[i]);
                local[9] = fmax(local[9], dataz[i]);
                local[10] += dataz[i];
                temp = pow(datax[i],2) + pow(datay[i],2) + pow(dataz[i],2);
                local[11] += temp;
                local[12] = fmax(temp, local[12]);
            }

            if(crank == 0)
            {
                piScalarData[0] = local[0];
                piScalarData[1] = local[1];
            }
            MPI_Reduce(local+2, piScalarData+2, 1, MPI_PRECISION, MPI_MIN, 0, ccomm);
            MPI_Reduce(local+3, piScalarData+3, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
            MPI_Reduce(local+4, piScalarData+4, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
            MPI_Reduce(local+5, piScalarData+5, 1, MPI_PRECISION, MPI_MIN, 0, ccomm);
            MPI_Reduce(local+6, piScalarData+6, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
            MPI_Reduce(local+7, piScalarData+7, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
            MPI_Reduce(local+8, piScalarData+8, 1, MPI_PRECISION, MPI_MIN, 0, ccomm);
            MPI_Reduce(local+9, piScalarData+9, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
            MPI_Reduce(local+10, piScalarData+10, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
            MPI_Reduce(local+11, piScalarData+11, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
            MPI_Reduce(local+12, piScalarData+12, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);

            if(magEquation)
            {
                datax = B->vec->x->spatial;
                datay = B->vec->y->spatial;
                dataz = B->vec->z->spatial;

                local[13] = datax[0];
                local[14] = datax[0];
                local[15] = datax[0];
                local[16] = datay[0];
                local[17] = datay[0];
                local[18] = datay[0];
                local[19] = dataz[0];
                local[20] = dataz[0];
                local[21] = dataz[0];
                temp = pow(datax[0],2) + pow(datay[0],2) + pow(dataz[0],2);
                local[22] = temp;
                local[23] = temp;
                for(i = 1; i < spatialCount; i++)
                {
                    local[13] = fmin(local[13], datax[i]);
                    local[14] = fmax(local[14], datax[i]);
                    local[15] += datax[i];
                    local[16] = fmin(local[16], datay[i]);
                    local[17] = fmax(local[17], datay[i]);
                    local[18] += datay[i];
                    local[19] = fmin(local[19], dataz[i]);
                    local[20] = fmax(local[20], dataz[i]);
                    local[21] += dataz[i];
                    temp = pow(datax[i],2) + pow(datay[i],2) + pow(dataz[i],2);
                    local[22] += temp;
                    local[23] = fmax(temp, local[21]);
                }

                MPI_Reduce(local+13, piScalarData+13, 1, MPI_PRECISION, MPI_MIN, 0, ccomm);
                MPI_Reduce(local+14, piScalarData+14, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
                MPI_Reduce(local+15, piScalarData+15, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
                MPI_Reduce(local+16, piScalarData+16, 1, MPI_PRECISION, MPI_MIN, 0, ccomm);
                MPI_Reduce(local+17, piScalarData+17, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
                MPI_Reduce(local+18, piScalarData+18, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
                MPI_Reduce(local+19, piScalarData+19, 1, MPI_PRECISION, MPI_MIN, 0, ccomm);
                MPI_Reduce(local+20, piScalarData+20, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
                MPI_Reduce(local+21, piScalarData+21, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
                MPI_Reduce(local+22, piScalarData+22, 1, MPI_PRECISION, MPI_SUM, 0, ccomm);
                MPI_Reduce(local+23, piScalarData+23, 1, MPI_PRECISION, MPI_MAX, 0, ccomm);
            }
            free(local);
            
            //Either store our data for later output, or perform the write to 
            //file
            if(crank == 0)
            {
                //make the means means and not sums
                piScalarData[4] /= (nx*ny*nz);
                piScalarData[7] /= (nx*ny*nz);
                piScalarData[10] /= (nx*ny*nz);
                if(magEquation)
                {
                    piScalarData[15] /= (nx*ny*nz);
                    piScalarData[18] /= (nx*ny*nz);
                    piScalarData[21] /= (nx*ny*nz);
                }

                scalarCount++;
                piScalarData += numScalar;

                //no reason to send this small amount of data to an IO node.
                //our root compute node will just quickly take care of this.
                if(scalarCount >= scalarPerF)
                {
                    char fileName[100];
                    sprintf(fileName, "Scalars/%08d",iteration);

                    FILE * out = fopen(fileName, "w");
                    fwrite(scalarData, sizeof(PRECISION), numScalar * scalarPerF, out);
                    fclose(out);

                    scalarCount = 0;
                    piScalarData = scalarData;
                }
            }
        }
    }
    
    //Time for spatial file output?
    if(iteration % spatialRate == 0)
    {
        char * name = (char *)malloc(100);
        
        //create the directory
        if(crank == 0)
        {
            sprintf(name, "Spatial/%08d",iteration);
            mkdir(name, S_IRWXU);

            //Record the simulation time that this snapshot belongs to
            sprintf(name, "Spatial/%08d/info",iteration);
            FILE * info;
            info = fopen(name, "w");
            fprintf(info, infostro, elapsedTime);
            fclose(info);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if(compute_node)
        {
            if(momEquation || kinematic)
            {
                trace("Outputing u\n");
                writeSpatial(u->vec->x,0);

                trace("Outputing v\n");
                writeSpatial(u->vec->y,0);

                trace("Outputing w\n");
                writeSpatial(u->vec->z,0);
            }
            if(tEquation)
            {
                trace("Outputting T\n");
                writeSpatial(T, 0);
            }
            if(magEquation)
            {
                trace("Outputing Bx\n");
                writeSpatial(B->vec->x,0);

                trace("Outputing By\n");
                writeSpatial(B->vec->y,0);

                trace("Outputing Bz\n");
                writeSpatial(B->vec->z,0);
            }
        }
        else if(io_node)
        {
            if(momEquation || kinematic)
            {
                sprintf(name, "Spatial/%08d/u", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);

                sprintf(name, "Spatial/%08d/v", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);

                sprintf(name, "Spatial/%08d/w", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);
            }
            if(tEquation)
            {
                sprintf(name, "Spatial/%08d/T", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);
            }
            if(magEquation)
            {
                sprintf(name, "Spatial/%08d/Bx", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);

                sprintf(name, "Spatial/%08d/By", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);

                sprintf(name, "Spatial/%08d/Bz", iteration);
                trace("Writing to file %s\n", name);
                writeSpatial(0, name);
            }
        }

        free(name);
    }

    if(iteration % checkRate == 0)
    {
        if(compute_node)
            writeCheckpoint();
    }
}


/*
 * This routine is virtually a memory dump.  There is no reason to gather data
 * to IO nodes, we are just going to open up files for each processor and dump
 * our section of the global array there.  We are also going to dump all of the
 * stored forcing evaluations.  The design is that any simulation that begins
 * from one of these checkpoints will be identical to a simulation that did
 * not stop in the first place.
 * 
 * Since this program is capable of terminating at any point in time, we take 
 * measures to ensure that the program does not terminate DURING a checkpoint
 * write and creating data corruption.  This is done by having two separate
 * checkpoint directories.  We alternate writes between these two directories,
 * and the final thing written during an output is a single file indicating 
 * that IO is done.  This way we always have a pristine checkpoint file that
 * we can restart from, regardless of when the program happens to terminate.
 * 
 * TODO: This folder gets very cluttered, especially with a lot of processors.
 *       Things should be changed so each processor gets one file, rather than
 *       one file per array being saved.
 */
void writeCheckpoint()
{
    FILE * out;
    char name[100];

    trace("Writing to Checkpoint%d\n", checkDir);

    if(momEquation)
    {
        sprintf(name,"Checkpoint%d/upol%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->poloidal->spectral, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/upolf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->poloidal->force1, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/upolf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->poloidal->force2, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/utor%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->toroidal->spectral, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/utorf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->toroidal->force1, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/utorf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->toroidal->force2, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umx%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_x, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umxf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_xf1, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umxf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_xf2, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umy%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_y, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umyf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_yf1, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umyf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_yf2, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umz%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(&(u->sol->mean_z), sizeof(complex<PRECISION>), 1, out);
        fclose(out);
    }

    if(magEquation)
    {
        sprintf(name,"Checkpoint%d/bpol%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->poloidal->spectral, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bpolf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->poloidal->force1, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bpolf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->poloidal->force2, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/btor%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->toroidal->spectral, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/btorf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->toroidal->force1, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/btorf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->toroidal->force2, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmx%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_x, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmxf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_xf1, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmxf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_xf2, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmy%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_y, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmyf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_yf1, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmyf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_yf2, sizeof(complex<PRECISION>), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmz%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(&(B->sol->mean_z), sizeof(complex<PRECISION>), 1, out);
        fclose(out);
    }

    if(tEquation)
    {
        sprintf(name,"Checkpoint%d/T%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(T->spectral, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/Tf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(T->force1, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/Tf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(T->force2, sizeof(complex<PRECISION>), spectralCount, out);
        fclose(out);
    }
    
    //finish all the important data.  Then update the state file so we know
    //things are completed
    MPI_Barrier(ccomm);

    if(crank == 0)
    {
        sprintf(name, "Checkpoint%d/state", checkDir);
        FILE * out = fopen(name, "w");
        fwrite(&elapsedTime, sizeof(PRECISION), 1, out);
        fwrite(&dt, sizeof(PRECISION), 1, out);
        fwrite(&dt1, sizeof(PRECISION), 1, out);
        fwrite(&iteration, sizeof(int), 1, out);
        fclose(out);
    }

    //Switch which directory we write to next
    if(checkDir == 0)
        checkDir = 1;
    else
        checkDir = 0;
}

/*
 * This is largely the inverse of the write method.  In order for this to work, 
 * the program MUST be run with the same number of compute nodes as previously
 * ran. 
 * 
 * The only extra bit is we here have to determine which checkpoint file is the
 * one to start from.  An earlier IO method ensures that the two directories 
 * exist and have a status file.  We first read in both status files, and the
 * the one with the latest simulation time is the newest and the one we will
 * begin from.  Since this status file is the last thing written and it is
 * written by a single processor after ALL of the other processors have finished
 * dumping data, we are guaranteed to not be reading in a checkpoint that was
 * corrupted because the program terminated while actually writing a checkpoint.
 */
void readCheckpoint()
{
    if(grank == 0)
    {
        PRECISION dElapsed;
        PRECISION dDt;
        PRECISION dDt1;
        int dIteration;

        FILE * out = fopen("Checkpoint0/state", "r");
        if(out == 0)
        {
            error("Failed to open Checpoint0 state file!  Crashing gracelessly...");
        }
        fread(&dElapsed, sizeof(PRECISION), 1, out);
        fread(&dDt, sizeof(PRECISION), 1, out);
        fread(&dDt1, sizeof(PRECISION), 1, out);
        fread(&dIteration, sizeof(int), 1, out);
        fclose(out);

        out = fopen("Checkpoint1/state", "r");
        if(out == 0)
        {
            error("Failed to open Checpoint0 state file!  Crashing gracelessly...");
        }
        fread(&elapsedTime, sizeof(PRECISION), 1, out);
        fread(&dt, sizeof(PRECISION), 1, out);
        fread(&dt1, sizeof(PRECISION), 1, out);
        fread(&iteration, sizeof(int), 1, out);
        fclose(out);

        if(iteration > dIteration)
        {
            checkDir = 1;
        }
        else
        {
            checkDir = 0;
            iteration = dIteration;
            elapsedTime = dElapsed;
            dt = dDt;
            dt1 = dDt1;
        }

    }

    //Let all processors know where we currently are in this simulation.
    MPI_Bcast(&iteration, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&elapsedTime, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt1, 1, MPI_PRECISION, 0, MPI_COMM_WORLD);
    MPI_Bcast(&checkDir, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);


    if(compute_node)
    {
        FILE * in;
        char name[100];

        trace("Reading from Checkpoint%d", checkDir);

        if(momEquation)
        {
            sprintf(name,"Checkpoint%d/upol%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->poloidal->spectral, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/upolf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->poloidal->force1, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/upolf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->poloidal->force2, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/utor%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->toroidal->spectral, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/utorf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->toroidal->force1, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/utorf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->toroidal->force2, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umx%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_x, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umxf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_xf1, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umxf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_xf2, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umy%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_y, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umyf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_yf1, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umyf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_yf2, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umz%d", checkDir, crank);
            in = fopen(name,"r");
            fread(&(u->sol->mean_z), sizeof(complex<PRECISION>), 1, in);
            fclose(in);

            recomposeSolenoidal(u->sol, u->vec);
            fftBackward(u->vec->x);
            fftBackward(u->vec->y);
            fftBackward(u->vec->z);
        }

        if(magEquation)
        {
            sprintf(name,"Checkpoint%d/bpol%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->poloidal->spectral, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bpolf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->poloidal->force1, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bpolf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->poloidal->force2, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/btor%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->toroidal->spectral, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/btorf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->toroidal->force1, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/btorf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->toroidal->force2, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmx%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_x, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmxf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_xf1, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmxf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_xf2, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmy%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_y, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmyf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_yf1, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmyf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_yf2, sizeof(complex<PRECISION>), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmz%d", checkDir, crank);
            in = fopen(name,"r");
            fread(&(B->sol->mean_z), sizeof(complex<PRECISION>), 1, in);
            fclose(in);

            recomposeSolenoidal(B->sol, B->vec);
            fftBackward(B->vec->x);
            fftBackward(B->vec->y);
            fftBackward(B->vec->z);
        }

        if(tEquation)
        {
            sprintf(name,"Checkpoint%d/T%d", checkDir, crank);
            in = fopen(name,"r");
            fread(T->spectral, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/Tf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(T->force1, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/Tf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(T->force2, sizeof(complex<PRECISION>), spectralCount, in);
            fclose(in);

            fftBackward(T);
        }

    }
}
