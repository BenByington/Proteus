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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include "IO.h"
#include "Field.h"
#include "Environment.h"
#include "Log.h"
#include "State.h"
#include "Numerics.h"
#include "Communication.h"

FILE * status = 0;

int scalarCount;
PRECISION * scalarData;
PRECISION * piScalarData;
int numScalar;

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
    writeSpatial(f, "Test/x");
    writeSpatial(f+1, "Test/y");
    writeSpatial(f+2, "Test/z");

    debug("Reading test IO data from files\n");
    readSpatial(f2, "Test/x");
    readSpatial(f2+1, "Test/y");
    readSpatial(f2+2, "Test/z");

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

void writeSpatial(field * f, char * name)
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

        displs[0] = 0;
        displs[1] = 0;
        rcvcounts[0] = 0;

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
    int disp = 0;
    for(i = 0; i < my_io_layer; i++)
    {
        for(j = io_layers[i].min; j <= io_layers[i].max; j++)
        {
            disp += all_z[j].width;
        }
    }
    disp *= nx * ny * sizeof(PRECISION);
    trace("Our view starts at element %d\n", disp);
    trace("Setting view...\n");
    MPI_File_set_view(fh, disp, MPI_PRECISION, MPI_PRECISION, "native", MPI_INFO_NULL);
    trace("Writing to file...\n");
    MPI_File_write(fh, sndbuff, nx * ny * nz_layers, MPI_PRECISION, MPI_STATUS_IGNORE );
    MPI_File_close(&fh);

    free(sndbuff);
    free(rcvbuff);

    debug("Write Spatial completed\n");

}

void readSpatial(field * f, char * name)
{
    int i,j,k,l,m;
    debug("Reading spatial data from file %s\n", name);

    PRECISION * sndbuff = 0;
    PRECISION * rcvbuff = 0;

    //check if there is an info file telling us our starting time
    if(grank == 0)
    {

    }

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
        MPI_File_set_view(fh, disp, MPI_PRECISION, MPI_PRECISION, "native", MPI_INFO_NULL);
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

        scalarData = malloc(numScalar * scalarPerF * sizeof(PRECISION));
        piScalarData = scalarData;
    }
}

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

#include "LogTrace.h"

void writeCheckpoint()
{
    FILE * out;
    char name[100];

    trace("Writing to Checkpoint%d\n", checkDir);

    if(momEquation)
    {
        sprintf(name,"Checkpoint%d/upol%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->poloidal->spectral, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/upolf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->poloidal->force1, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/upolf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->poloidal->force2, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/utor%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->toroidal->spectral, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/utorf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->toroidal->force1, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/utorf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->toroidal->force2, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umx%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_x, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umxf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_xf1, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umxf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_xf2, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umy%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_y, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umyf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_yf1, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umyf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(u->sol->mean_yf2, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/umz%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(&(u->sol->mean_z), sizeof(complex PRECISION), 1, out);
        fclose(out);
    }

    if(magEquation)
    {
        sprintf(name,"Checkpoint%d/bpol%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->poloidal->spectral, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bpolf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->poloidal->force1, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bpolf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->poloidal->force2, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/btor%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->toroidal->spectral, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/btorf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->toroidal->force1, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/btorf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->toroidal->force2, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmx%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_x, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmxf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_xf1, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmxf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_xf2, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmy%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_y, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmyf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_yf1, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmyf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(B->sol->mean_yf2, sizeof(complex PRECISION), ndkz, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/bmz%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(&(B->sol->mean_z), sizeof(complex PRECISION), 1, out);
        fclose(out);
    }

    if(tEquation)
    {
        sprintf(name,"Checkpoint%d/T%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(T->spectral, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/Tf1%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(T->force1, sizeof(complex PRECISION), spectralCount, out);
        fclose(out);

        sprintf(name,"Checkpoint%d/Tf2%d", checkDir, crank);
        out = fopen(name,"w");
        fwrite(T->force2, sizeof(complex PRECISION), spectralCount, out);
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
            fread(u->sol->poloidal->spectral, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/upolf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->poloidal->force1, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/upolf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->poloidal->force2, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/utor%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->toroidal->spectral, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/utorf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->toroidal->force1, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/utorf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->toroidal->force2, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umx%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_x, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umxf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_xf1, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umxf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_xf2, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umy%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_y, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umyf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_yf1, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umyf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(u->sol->mean_yf2, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/umz%d", checkDir, crank);
            in = fopen(name,"r");
            fread(&(u->sol->mean_z), sizeof(complex PRECISION), 1, in);
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
            fread(B->sol->poloidal->spectral, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bpolf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->poloidal->force1, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bpolf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->poloidal->force2, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/btor%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->toroidal->spectral, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/btorf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->toroidal->force1, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/btorf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->toroidal->force2, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmx%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_x, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmxf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_xf1, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmxf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_xf2, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmy%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_y, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmyf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_yf1, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmyf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(B->sol->mean_yf2, sizeof(complex PRECISION), ndkz, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/bmz%d", checkDir, crank);
            in = fopen(name,"r");
            fread(&(B->sol->mean_z), sizeof(complex PRECISION), 1, in);
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
            fread(T->spectral, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/Tf1%d", checkDir, crank);
            in = fopen(name,"r");
            fread(T->force1, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            sprintf(name,"Checkpoint%d/Tf2%d", checkDir, crank);
            in = fopen(name,"r");
            fread(T->force2, sizeof(complex PRECISION), spectralCount, in);
            fclose(in);

            fftBackward(T);
        }

    }
}
