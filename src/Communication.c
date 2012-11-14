/* 
 * File:   Communication.cpp
 * Author: Ben
 * 
 * Created on February 18, 2010, 1:15 PM
 */

#include "Communication.h"
#include "mpi.h"

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "Environment.h"
#include "Log.h"

int whichfft;

void clean();

fftw_plan planf1;
fftw_plan planf2;
fftw_plan planf3;
fftw_plan planb1;
fftw_plan planb2;
fftw_plan planb3;

void initfft1();
void fft1_forward(double * in, complex double * out);
void fft1_tpf1(complex double * in, complex double * out);
void fft1_tpf2(complex double * in, complex double * out);
void fft1_backward(complex double * in, double * out);
void fft1_tpb1(complex double * in, complex double * out);
void fft1_tpb2(complex double * in, complex double * out);

void fft_tpf3(complex double * in, complex double * out);
void fft_tpb3(complex double * in, complex double * out);

void initfft2();
void fft2_forward(double * in, complex double * out);
void fft2_tpf1(complex double * in, complex double * out);
void fft2_tpf2(complex double * in, complex double * out);
void fft2_backward(complex double * in, double * out);
void fft2_tpb1(complex double * in, complex double * out);
void fft2_tpb2(complex double * in, complex double * out);

void testfft1();
void repeatfft1(int count);
void testfft2();
void repeatfft2(int count);

void generateFunc(int * ks, int len, double * out);

void com_init(int measure)
{
    int test = 1;

    info("Initializing FFT routines.  Measure = %d\n", measure);
    if(measure)
    {

        clock_t start;
        clock_t stop;

        initfft1();
        if(test)
        {   
            info("Testing fft1\n",0);
            testfft1();
        }

        debug("Timing fft1\n",0);
        start = clock();
        repeatfft1(100);
        stop = clock();
        int diff1 = stop - start;

        fftw_cleanup();
        initfft2();
        if(test)
        {
            info("Testing fft2\n",0);
            testfft2();
        }

        debug("Timing fft2\n",0);
        start = clock();
        repeatfft2(100);
        stop = clock();

        int diff2 = stop - start;

        if(test)
        {
            info("fft1 finished in %g\n", (double)diff1 / CLOCKS_PER_SEC);
            info("fft2 finished in %g\n", (double)diff2 / CLOCKS_PER_SEC);
        }

        if(diff1 < diff2)
        {
            fftw_cleanup();
            initfft1();
            whichfft = FFT1;
        }
        else
        {
            whichfft = FFT2;
        }
    }
    else
    {
        initfft1();
        whichfft = FFT1;
    }

    info("FFT %d is in use for this run\n", whichfft);
}

/*
 * This inits the code for the 3Dfft that takes data in the [z][x][y] layout.
 * All array packing and transposing is handled explicitly, and fftw only does
 * contiguous transforms
 */
void initfft1()
{
    debug("Initializing fft1...\n",0);
    double * real;
    complex double * comp1;
    complex double * comp2;

    real = (double*)fftw_malloc(my_z->width * my_x->width * ny * sizeof(double));
    comp1 = (complex double*)fftw_malloc(my_z->width * my_x->width * nky * sizeof(complex double));
    planf1 = fftw_plan_many_dft_r2c(1, &ny, my_x->width * my_z->width, (double *)real, 0, 1, ny, comp1, 0, 1, nky, FFTW_MEASURE);
    planb1 = fftw_plan_many_dft_c2r(1, &ny, my_x->width * my_z->width, comp1, 0, 1, nky, (double*)real, 0, 1, ny, FFTW_MEASURE);
    fftw_free(real);
    fftw_free(comp1);

    comp1 = (complex double*)fftw_malloc(my_z->width * my_ky->width * nx * sizeof(complex double));
    comp2 = (complex double*)fftw_malloc(my_z->width * my_ky->width * nkx * sizeof(complex double));
    planf2 = fftw_plan_many_dft(1, &nx, my_z->width * my_ky->width, comp1, 0, 1, nx, comp2, 0, 1, nkx, FFTW_FORWARD, FFTW_MEASURE);
    planb2 = fftw_plan_many_dft(1, &nx, my_z->width * my_ky->width, comp2, 0, 1, nkx, comp1, 0, 1, nx, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_free(comp1);
    fftw_free(comp2);

    comp1 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nz * sizeof(complex double));
    comp2 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex double));
    planf3 = fftw_plan_many_dft(1, &nz, my_kx->width * my_ky->width, comp1, 0, 1, nz, comp2, 0, 1, nkz, FFTW_FORWARD, FFTW_MEASURE);
    planb3 = fftw_plan_many_dft(1, &nz, my_kx->width * my_ky->width, comp2, 0, 1, nkz, comp1, 0, 1, nz, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_free(comp1);
    fftw_free(comp2);

    debug("Initialization done\n",0);
}

void fft1_forward(double * in, complex double* out)
{
    int i;

    trace("Begin fftw1 forward transform\n",0);
    int mySize1 = my_z->width * my_x->width * nky;
    int mySize2 = my_z->width * my_ky->width * nx;
    int mySize3 = my_kx->width * my_ky->width * nz;

    complex double * comp = (complex double*)fftw_malloc(mySize1*sizeof(complex double));
    fftw_execute_dft_r2c(planf1, in, comp);

    complex double * comp2 = (complex double*)fftw_malloc(mySize2*sizeof(complex double));
    fft1_tpf1(comp, comp2);
    fftw_free(comp);

    complex double * comp3 = (complex double*)fftw_malloc(mySize2*sizeof(complex double));
    fftw_execute_dft(planf2, comp2, comp3);
    fftw_free(comp2);

    complex double * comp4 = (complex double*)fftw_malloc(mySize3*sizeof(complex double));
    fft1_tpf2(comp3, comp4);
    fftw_free(comp3);

    complex double * comp5 = (complex double*)fftw_malloc(mySize3*sizeof(complex double));
    fftw_execute_dft(planf3, comp4, comp5);
    fftw_free(comp4);

    fft_tpf3(comp5, out);
    fftw_free(comp5);

    int size = my_kx->width * my_ky->width * ndkz;
    double factor = sqrt(ny) * sqrt(nx) * sqrt(nz);
    for(i = 0; i < size; i++)
        out[i] /= factor;

    trace("Forward fftw competed\n",0);
}

void fft1_tpf1(complex double* in, complex double* out)
{
    int i,j,k;
    int maxSize1 = max_z->width * max_x->width * max_ky->width;
    trace("Starting first transpose for forward fft1\n",0);
    //in is complex double[my_z->width][my_x->width][nky]
    //sndbuff and rcvbuff is complex double[hsize][maxz->width][maxx->width][maxky->width]
    complex double * sndbuff = (complex double*)malloc(maxSize1 * hsize * sizeof(complex double));
    complex double * rcvbuff = (complex double*)malloc(maxSize1 * hsize * sizeof(complex double));

    complex double * piin = in;
    complex double * pisbuff = sndbuff;

    trace("Packing arrays for MPI all-to-all\n",0);
    //loop over each x and z to process contiguous 1D arrays
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_x->width; j++)
        {
            //pisbuff isn't moved contiguously through memory.  It needs
            //to be set here for the processing of each array
            pisbuff = sndbuff + j * max_ky->width + i * max_ky->width * max_x->width;

            //loop over the processors we need to divide this array among
            for(k = 0; k < hsize; k++)
            {
                //our internal pointers start where they need to.  Just do the
                //memcpy
                memcpy(pisbuff, piin, all_ky[k].width * sizeof(complex double));

                //incriment pointers.  Simple for picomp.  Be careful for
                //pisbuff
                piin += all_ky[k].width;
                pisbuff += maxSize1;
            }
            //There is a tail of wavelengths we are skipping over at the end
            //for dealiasing.  Adjust the internal pointer to reflect this
            piin += dealias_ky.width;
        }
    }

    trace("Sending data over network\n",0);
    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_DOUBLE, rcvbuff, 2 * maxSize1, MPI_DOUBLE, hcomm);

    trace("Unpacking data from transfer\n",0);
    //out is complex double[myz->width][my_ky->width][nx]
    //rcvbuff is complex double[hsize][maxz->width][maxx->width][maxky->width]
    //loop overy every element in out and set it to the correct value.  This will
    //perform the transpost
    //TODO look into a more efficient way to do this...
    int xBig;
    int xSmall;
    complex double * pirbuff = rcvbuff;
    complex double * piout = out;
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            xBig = 0;
            for(k = 0; k < nx; k++)
            {
                //find the current slot for pirbuff.
                if(k > all_x[xBig].max)
                    xBig++;
                xSmall = k - all_x[xBig].min;
                pirbuff = rcvbuff + j + xSmall * max_ky->width + i * max_x->width * max_ky->width + xBig * max_z->width * max_x->width * max_ky->width;

                (*piout) = (*pirbuff);

                //move piout to the next slot
                piout++;
            }
        }
    }

    free(sndbuff);
    free(rcvbuff);
}

void fft1_tpf2(complex double* in, complex double* out)
{
    trace("Starting second transpose for forward fft1\n",0);
    int i,j,k;
    int maxSize1 = max_z->width * max_kx->width * max_ky->width;

    //in is complex double[my_z->width][my_ky->width][nkx]
    //sndbuff and rcvbuff is complex double[vsize][max_z->width][max_ky->width][max_kx->width]
    complex double * sndbuff = (complex double*)malloc(maxSize1 * vsize * sizeof(complex double));
    complex double * rcvbuff = (complex double*)malloc(maxSize1 * vsize * sizeof(complex double));

    complex double * piin = in;
    complex double * pisbuff = sndbuff;

    //One of the processors is going to have to deal with skipping over wavelengths
    //for dealiasing.  Stay tuned to find out who!!
    int dProc = -1;
    int nlow = -1;
    int nhigh = -1;
    for(i = 0; i < vsize; i++)
    {
        if(all_kx[i].min <= dealias_kx.min && all_kx[i].max >= dealias_kx.min)
        {
            dProc = i;
            nlow = dealias_kx.min - all_kx[i].min;
            nhigh = all_kx[i].width - nlow;
        }
    }

    trace("Packing arrays for MPI all-to-all\n",0);
    //loop over each y and z to process contiguous 1D arrays
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            //pisbuff isn't moved contiguously through memory.  It needs
            //to be set here for the processing of each array
            pisbuff = sndbuff + j * max_kx->width + i * max_ky->width * max_kx->width;

            //loop over the processors we need to divide this array among
            for(k = 0; k < vsize; k++)
            {
                //handle the skipping of dealiase wavelengths
                if(k == dProc)
                {
                    
                    if(nlow)
                        memcpy(pisbuff, piin, nlow * sizeof(complex double));

                    piin += nlow + dealias_kx.width;

                    if(nhigh)
                        memcpy(pisbuff + nlow, piin, nhigh * sizeof(complex double));
                    piin += nhigh;
                }
                else
                {
                    //our internal pointers start where they need to.  Just do the
                    //memcpy
                    memcpy(pisbuff, piin, all_kx[k].width * sizeof(complex double));

                    //incriment pointers.  Simple for picomp.  Be careful for
                    //pisbuff
                    piin += all_kx[k].width;
                }
                pisbuff += maxSize1;
            }
        }
    }

    trace("Sending data over network\n",0);
    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_DOUBLE, rcvbuff, 2 * maxSize1, MPI_DOUBLE, vcomm);

    trace("Unpacking data from transfer\n",0);
    //out is complex double[my_kx->width][my_ky->width][nz]
    //rcvbuff is complex double[hsize][max_z->width][max_ky->width][max_kx->width]
    //loop overy every element in out and set it to the correct value.  This will
    //perform the transpost
    //TODO look into a more efficient way to do this...
    int zBig;
    int zSmall;
    complex double * pirbuff = rcvbuff;
    complex double * piout = out;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            zBig = 0;
            for(k = 0; k < nz; k++)
            {
                //find the current slot for pirbuff.
                if(k > all_z[zBig].max)
                    zBig++;
                zSmall = k - all_z[zBig].min;
                pirbuff = rcvbuff + i + j * max_kx->width + zSmall * max_ky->width * max_kx->width + zBig * max_z->width * max_ky->width * max_kx->width;

                (*piout) = (*pirbuff);

                //move piout to the next slot
                piout++;
            }
        }
    }

    free(sndbuff);
    free(rcvbuff);
}

void fft_tpf3(complex double * in, complex double * out)
{
    trace("Performing final dealias for forward transform\n",0);
    int i,j,k;
    //in is [my_kx][my_ky][nkz]
    //out is [my_kx][my_ky][ndkz]

    complex double * piin = in;
    complex double * piout = out;
    //loop over x and y to process contiguous arrays

    int len1 = dealias_kz.min;
    int cut = dealias_kz.width;
    int len2 = nkz - len1 - cut;

    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            memcpy(piout, piin, len1 * sizeof(complex double));

            piout += len1;
            piin += len1 + cut;

            memcpy(piout, piin, len2*sizeof(complex double));

            piout += len2;
            piin += len2;
        }
    }
}

void fft1_backward(complex double* in, double * out)
{
    trace("Begin fft1 backwards transform\n",0);
    int mySize1 = my_kx->width * my_ky->width * nz;
    int mySize2 = my_z->width * my_ky->width * nkx;
    int mySize3 = my_z->width * my_x->width * nky;

    complex double * comp = (complex double*)fftw_malloc(mySize1 * sizeof(complex double));
    fft_tpb3(in, comp);

    complex double * comp2 = (complex double*)fftw_malloc(mySize1 * sizeof(complex double));
    fftw_execute_dft(planb3, comp, comp2);
    fftw_free(comp);

    complex double * comp3 = (complex double*)fftw_malloc(mySize2 * sizeof(complex double));
    fft1_tpb2(comp2, comp3);
    fftw_free(comp2);

    complex double * comp4 = (complex double*)fftw_malloc(mySize2 * sizeof(complex double));
    fftw_execute_dft(planb2, comp3, comp4);
    fftw_free(comp3);

    complex double * comp5 = (complex double*)fftw_malloc(mySize3 * sizeof(complex double));
    fft1_tpb1(comp4, comp5);
    fftw_free(comp4);

    fftw_execute_dft_c2r(planb1, comp5, out);
    fftw_free(comp5);
    

    int i;
    double factor = sqrt(ny) * sqrt(nx) * sqrt(nz);
    int size = my_x->width * my_z->width * ny;
    for(i = 0; i < size; i++)
    {
        out[i] /= factor;
    }

    trace("Inverse FFT completed\n",0);
}

void fft1_tpb1(complex double* in, complex double* out)
{
    int i,j,k;
    //in is complex double[my_z->width][my_ky->width][nx]
    //sndbuff is complex double[hsize][maxz->width][maxx->width][maxky->width]
    int maxSize1 = max_z->width * max_x->width * max_ky->width;
    complex double* sndbuff = (complex double*)malloc(maxSize1 * hsize * sizeof(complex double));
    complex double* rcvbuff = (complex double*)malloc(maxSize1 * hsize * sizeof(complex double));


    //loop overy every element in in and send it to the correct spot in sndbuff.
    //This will perform the transpose
    //TODO look into a more efficient way to do this...
    int xBig;
    int xSmall;
    complex double * pisbuff = sndbuff;
    complex double * piin = in;
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            xBig = 0;
            for(k = 0; k < nx; k++)
            {
                //find the current slot for pisbuff.
                if(k > all_x[xBig].max)
                    xBig++;
                xSmall = k - all_x[xBig].min;
                pisbuff = sndbuff + j + xSmall * max_ky->width + i * max_x->width * max_ky->width + xBig * max_z->width * max_x->width * max_ky->width;

                (*pisbuff) = (*piin);

                //move piout to the next slot
                piin++;
            }
        }
    }

    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_DOUBLE, rcvbuff, 2 * maxSize1, MPI_DOUBLE, hcomm);

    complex double * pirbuff = rcvbuff;
    complex double * piout = out;
    //loop over each x and z to process contiguous 1D arrays
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_x->width; j++)
        {
            //pisbuff isn't moved contiguously through memory.  It needs
            //to be set here for the processing of each array
            pirbuff = rcvbuff + j * max_ky->width + i * max_ky->width * max_x->width;

            //loop over the processors
            for(k = 0; k < hsize; k++)
            {
                //our internal pointers start where they need to.  Just do the
                //memcpy
                memcpy(piout, pirbuff, all_ky[k].width * sizeof(complex double));

                piout += all_ky[k].width;
                pirbuff += maxSize1;
            }
            //There are some dealiased wavelengths at the end that we need
            //to put back in as 0's
            memset(piout, 0, dealias_ky.width * sizeof(complex double));
            piout += dealias_ky.width;
        }
    }

    free(sndbuff);
    free(rcvbuff);
}

void fft1_tpb2(complex double* in, complex double* out)
{
    int i,j,k;
    int maxSize1 = max_z->width * max_kx->width * max_ky->width;
    complex double * sndbuff = (complex double*)malloc(maxSize1 * vsize * sizeof(complex double));
    complex double * rcvbuff = (complex double*)malloc(maxSize1 * vsize * sizeof(complex double));


    //in is complex double[my_kx->width][my_ky->width][nz]
    //sndbuff is complex double[vsize][max_z->width][max_ky->width][max_kx->width]
    //loop overy every element in in and set it to the correct value.  This will
    //perform the transpost
    //TODO look into a more efficient way to do this...
    int zBig;
    int zSmall;
    complex double * pisbuff = sndbuff;
    complex double * piin = in;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            zBig = 0;
            for(k = 0; k < nz; k++)
            {
                //find the current slot for pirbuff.
                if(k > all_z[zBig].max)
                    zBig++;
                zSmall = k - all_z[zBig].min;
                pisbuff = sndbuff + i + j * max_kx->width + zSmall * max_ky->width * max_kx->width + zBig * max_z->width * max_ky->width * max_kx->width;

                (*pisbuff) = (*piin);

                //move piout to the next slot
                piin++;
            }
        }
    }

    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_DOUBLE, rcvbuff, 2 * maxSize1, MPI_DOUBLE, vcomm);

    complex double * piout = out;
    complex double * pirbuff = sndbuff;

    //One of the processors is going to have to deal with skipping over wavelengths
    //for dealiasing.  Stay tuned to find out who!!
    int dProc = -1;
    int nlow = -1;
    int nhigh = -1;
    for(i = 0; i < vsize; i++)
    {
        if(all_kx[i].min <= dealias_kx.min && all_kx[i].max >= dealias_kx.min)
        {
            dProc = i;
            nlow = dealias_kx.min - all_kx[i].min;
            nhigh = all_kx[i].width - nlow;
        }
    }

    //loop over each y and z to process contiguous 1D arrays
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            //pisbuff isn't moved contiguously through memory.  It needs
            //to be set here for the processing of each array
            pirbuff = rcvbuff + j * max_kx->width + i * max_ky->width * max_kx->width;

            //loop over the processors we need to divide this array among
            for(k = 0; k < vsize; k++)
            {
                //handle the skipping of dealiase wavelengths
                if(k == dProc)
                {

                    if(nlow)
                        memcpy(piout, pirbuff, nlow * sizeof(complex double));

                    piout += nlow;
                    memset(piout, 0, dealias_kx.width*sizeof(complex double));
                    piout += dealias_kx.width;

                    if(nhigh)
                        memcpy(piout, pirbuff + nlow, nhigh * sizeof(complex double));
                    piout += nhigh;
                }
                else
                {
                    //our internal pointers start where they need to.  Just do the
                    //memcpy
                    memcpy(piout, pirbuff, all_kx[k].width * sizeof(complex double));

                    //incriment pointers.  Simple for picomp.  Be careful for
                    //pisbuff
                    piout += all_kx[k].width;
                }
                pirbuff += maxSize1;
            }
        }
    }

    free(sndbuff);
    free(rcvbuff);
}

void fft_tpb3(complex double * in, complex double * out)
{
    int i,j;
    //in is [my_kx][my_ky][ndkz]
    //out is [my_kx][my_ky][nkz]

    complex double * piin = in;
    complex double * piout = out;
    //loop over x and y to process contiguous arrays

    int len1 = dealias_kz.min;
    int cut = dealias_kz.width;
    int len2 = nkz - len1 - cut;

    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            memcpy(piout, piin, len1 * sizeof(complex double));
            piout += len1;
            
            memset(piout, 0, cut * sizeof(complex double));
            piout += cut;

            piin += len1;
            memcpy(piout, piin, len2*sizeof(complex double));

            piin += len2;
            piout += len2;
        }
    }
}

void initfft2()
{
    double * real;
    complex double * comp1;
    complex double * comp2;

    real = (double*)fftw_malloc(my_z->width * my_x->width * ny * sizeof(double));
    comp1 = (complex double*)fftw_malloc(nky * my_z->width * my_x->width * sizeof(complex double));
    planf1 = fftw_plan_many_dft_r2c(1, &ny, my_x->width * my_z->width, (double *)real, 0, 1, ny, comp1, 0, my_x->width * my_z->width, 1, FFTW_MEASURE);
    planb1 = fftw_plan_many_dft_c2r(1, &ny, my_x->width * my_z->width, comp1, 0, my_x->width * my_z->width, 1, (double*)real, 0, 1, ny, FFTW_MEASURE);
    fftw_free(real);
    fftw_free(comp1);

    comp1 = (complex double*)fftw_malloc(my_z->width * my_ky->width * nx * sizeof(complex double));
    comp2 = (complex double*)fftw_malloc(my_z->width * my_ky->width * nkx * sizeof(complex double));
    planf2 = fftw_plan_many_dft(1, &nx, my_z->width * my_ky->width, comp1, 0, 1, nx, comp2, 0, my_z->width * my_ky->width, 1, FFTW_FORWARD, FFTW_MEASURE);
    planb2 = fftw_plan_many_dft(1, &nx, my_z->width * my_ky->width, comp2, 0, my_z->width * my_ky->width, 1, comp1, 0, 1, nx, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_free(comp1);
    fftw_free(comp2);

    comp1 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nz * sizeof(complex double));
    comp2 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex double));
    planf3 = fftw_plan_many_dft(1, &nz, my_kx->width * my_ky->width, comp1, 0, 1, nz, comp2, 0, 1, nkz, FFTW_FORWARD, FFTW_MEASURE);
    planb3 = fftw_plan_many_dft(1, &nz, my_kx->width * my_ky->width, comp2, 0, 1, nkz, comp1, 0, 1, nz, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_free(comp1);
    fftw_free(comp2);
}

void fft2_forward(double* in, complex double* out)
{
    complex double * comp1 = (complex double*)fftw_malloc(nky * my_z->width * my_x->width * sizeof(complex double));
    fftw_execute_dft_r2c(planf1, in, comp1);

    complex double * comp2 = (complex double*)fftw_malloc(my_ky->width * my_z->width * nx * sizeof(complex double));
    fft2_tpf1(comp1, comp2);
    fftw_free(comp1);

    complex double * comp3 = (complex double*)fftw_malloc(nkx * my_ky->width * my_z->width * sizeof(complex double));
    fftw_execute_dft(planf2, comp2, comp3);
    fftw_free(comp2);

    complex double* comp4 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nz * sizeof(complex double));
    fft2_tpf2(comp3, comp4);
    fftw_free(comp3);

    complex double * comp5 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex double));
    fftw_execute_dft(planf3, comp4, comp5);
    fftw_free(comp4);

    fft_tpf3(comp5, out);
    fftw_free(comp5);

    int i;
    double factor = sqrt(ny) * sqrt(nx) * sqrt(nz);
    int size = my_kx->width * my_ky->width * ndkz;
    for(i = 0; i < size; i++)
        out[i] /= factor;
}

void fft2_tpf1(complex double* in, complex double* out)
{
    int i,j,k;
    complex double * rcvbuff = (complex double*)malloc(my_ky->width * my_z->width * nx * sizeof(complex double));
    //set up the send/receive data structures
    int * scnt = (int*)malloc(hsize * sizeof(int));
    int * sdisp = (int*)malloc(hsize * sizeof(int));
    int * rcnt = (int*)malloc(hsize * sizeof(int));
    int * rdisp = (int*)malloc(hsize * sizeof(int));

    scnt[0] = 2 * all_ky[0].width * my_z->width * my_x->width;
    rcnt[0] = 2 * my_ky->width * my_z->width * all_x[0].width;
    sdisp[0] = 0;
    rdisp[0] = 0;
    for(i = 1; i < hsize; i++)
    {
        scnt[i] = 2 * all_ky[i].width * my_z->width * my_x->width;
        rcnt[i] = 2 * my_ky->width * my_z->width * all_x[i].width;
        sdisp[i] = sdisp[i-1] + scnt[i-1];
        rdisp[i] = rdisp[i-1] + rcnt[i-1];
    }

    MPI_Alltoallv(in, scnt, sdisp, MPI_DOUBLE, rcvbuff, rcnt, rdisp, MPI_DOUBLE, hcomm);

    //rcvbuff has a very non-uniform layout, so we will simply things by moving
    //contiguously through it, and jumping around in out.
    //out = [my_ky->width][my_z->width][nx]
    //rcvbuff = [p][my_ky->width][my_z->width][px]
    complex double * pirbuff = rcvbuff;
    complex double *  piout = out;

    int offset = 0;
    for(i = 0; i < hsize; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < my_z->width; k++)
            {
                piout = out + offset + k * nx + j * my_z->width * nx;
                memcpy(piout, pirbuff, all_x[i].width * sizeof(complex double));
                pirbuff += all_x[i].width;
            }
        }
        offset += all_x[i].width;
    }

    free(rcvbuff);
    free(rcnt);
    free(rdisp);
    free(scnt);
    free(sdisp);
}

void fft2_tpf2(complex double* in, complex double* out)
{
    int i,j,k;
    complex double * rcvbuff = (complex double*)malloc(my_kx->width * my_ky->width * nz * sizeof(complex double));
    //set up the send/receive data structures
    int * scnt = (int*)malloc(vsize * sizeof(int));
    int * sdisp = (int*)malloc(vsize * sizeof(int));
    int * rcnt = (int*)malloc(vsize * sizeof(int));
    int * rdisp = (int*)malloc(vsize * sizeof(int));

    //For one of the processors, the information needed is not contiguous because
    //it is interrupted by wavelengths we wish to discard for dealiasing.
    int dProc = -1;
    int nlow = -1;
    int nhigh = -1;
    for(i = 0; i < vsize; i++)
    {
        if(all_kx[i].min <= dealias_kx.min && all_kx[i].max >= dealias_kx.min)
        {
            dProc = i;
            nlow = dealias_kx.min - all_kx[i].min;
            nhigh = all_kx[i].width - nlow;

            //fprintf(stderr, "Proc %d: %d %d %d %d %d %d  %d\n", i, all_kx[i].min, all_kx[i].max, dealias_kx.min, dealias_kx.max, dealias_kx.width, nlow, nhigh);
        }
    }
    //now move the data so the data for dProc is contiguous
    memmove(in + dealias_kx.min * my_ky->width * my_z->width, in + (dealias_kx.max+1) * my_ky->width * my_z->width, nhigh * my_ky->width * my_z->width * sizeof(complex double));

    scnt[0] = 2 * all_kx[0].width * my_ky->width * my_z->width;
    rcnt[0] = 2 * my_kx->width * my_ky->width * all_z[0].width;
    sdisp[0] = 0;
    rdisp[0] = 0;
    for(i = 1; i < vsize; i++)
    {
        scnt[i] = 2 * all_kx[i].width * my_ky->width * my_z->width;
        rcnt[i] = 2 * my_kx->width * my_ky->width * all_z[i].width;
        sdisp[i] = sdisp[i-1] + scnt[i-1];
        rdisp[i] = rdisp[i-1] + rcnt[i-1];

        //If we are the proc after dProc, we actually have a bunch of empty
        //wave modes to skip over...
        if(i == dProc+1)
        {
            sdisp[i] += 2 * dealias_kx.width * my_ky->width * my_z->width;
        }
    }

    MPI_Alltoallv(in, scnt, sdisp, MPI_DOUBLE, rcvbuff, rcnt, rdisp, MPI_DOUBLE, vcomm);

    //rcvbuff has a very non-uniform layout, so we will simplify things by moving
    //contiguously through it, and jumping around in out.
    //out = [my_kx->width][my_ky->width][nz]
    //rcvbuff = [p][my_kx->width][my_ky->width][pz]
    complex double * pirbuff = rcvbuff;
    complex double *  piout = out;

    int offset = 0;
    for(i = 0; i < vsize; i++)
    {
        for(j = 0; j < my_kx->width; j++)
        {
            for(k = 0; k < my_ky->width; k++)
            {
                piout = out + offset + k * nz + j * my_ky->width * nz;
                memcpy(piout, pirbuff, all_z[i].width * sizeof(complex double));
                pirbuff += all_z[i].width;
            }
        }
        offset += all_z[i].width;
    }

    free(rcvbuff);
    free(rcnt);
    free(rdisp);
    free(scnt);
    free(sdisp);
}

void fft2_backward(complex double* in, double* out)
{
    complex double * comp = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex double));
    fft_tpb3(in, comp);

    complex double * comp1 = (complex double*)fftw_malloc(my_kx->width * my_ky->width * nkz*sizeof(complex double));
    fftw_execute_dft(planb3, comp, comp1);
    fftw_free(comp);

    complex double * comp2 = (complex double*)fftw_malloc(nkx * my_ky->width * my_z->width*sizeof(complex double));
    fft2_tpb2(comp1, comp2);
    fftw_free(comp1);

    complex double * comp3 = (complex double*)fftw_malloc(my_ky->width * my_z->width * nx*sizeof(complex double));
    fftw_execute_dft(planb2, comp2, comp3);
    fftw_free(comp2);
    
    complex double * comp4 = (complex double*)fftw_malloc(nky * my_z->width * my_x->width*sizeof(complex double));
    fft2_tpb1(comp3, comp4);
    fftw_free(comp3);
    
    fftw_execute_dft_c2r(planb1, comp4, out);
    fftw_free(comp4);

    int i;
    double factor = sqrt(ny) * sqrt(nx) * sqrt(nz);
    int size = my_x->width * my_z->width * ny;
    for(i = 0; i < size; i++)
        out[i] /= factor;
}

void fft2_tpb1(complex double* in, complex double* out)
{
    int i,j,k;
    complex double * sndbuff = (complex double*)malloc(my_ky->width * my_z->width * nx * sizeof(complex double));
    //set up the send/receive data structures
    int * scnt = (int*)malloc(hsize * sizeof(int));
    int * sdisp = (int*)malloc(hsize * sizeof(int));
    int * rcnt = (int*)malloc(hsize * sizeof(int));
    int * rdisp = (int*)malloc(hsize * sizeof(int));

    //sndbuff has a very non-uniform layout, so we will simply things by moving
    //contiguously through it, and jumping around in out.
    //out = [my_ky->width][my_z->width][nx]
    //rcvbuff = [p][my_ky->width][my_z->width][px]
    complex double * pisbuff = sndbuff;
    complex double *  piin = in;

    int offset = 0;
    for(i = 0; i < hsize; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < my_z->width; k++)
            {
                piin = in + offset + k * nx + j * my_z->width * nx;
                memcpy(pisbuff, piin, all_x[i].width * sizeof(complex double));
                pisbuff += all_x[i].width;
            }
        }
        offset += all_x[i].width;
    }

    rcnt[0] = 2 * all_ky[0].width * my_z->width * my_x->width;
    scnt[0] = 2 * my_ky->width * my_z->width * all_x[0].width;
    sdisp[0] = 0;
    rdisp[0] = 0;
    for(i = 1; i < hsize; i++)
    {
        rcnt[i] = 2 * all_ky[i].width * my_z->width * my_x->width;
        scnt[i] = 2 * my_ky->width * my_z->width * all_x[i].width;
        sdisp[i] = sdisp[i-1] + scnt[i-1];
        rdisp[i] = rdisp[i-1] + rcnt[i-1];
    }

    MPI_Alltoallv(sndbuff, scnt, sdisp, MPI_DOUBLE, out, rcnt, rdisp, MPI_DOUBLE, hcomm);

    //make sure the dealiased wavelengths are 0
    memset(out + dealias_ky.min * my_x->width * my_z->width, 0, dealias_ky.width * my_x->width * my_z->width * sizeof(complex double));

    free(sndbuff);
    free(rcnt);
    free(rdisp);
    free(scnt);
    free(sdisp);
}

void fft2_tpb2(complex double* in, complex double* out)
{
    int i,j,k;
    complex double * sndbuff = (complex double*)malloc(my_kx->width * my_ky->width * nz * sizeof(complex double));
    //set up the send/receive data structures
    int * scnt = (int*)malloc(vsize * sizeof(int));
    int * sdisp = (int*)malloc(vsize * sizeof(int));
    int * rcnt = (int*)malloc(vsize * sizeof(int));
    int * rdisp = (int*)malloc(vsize * sizeof(int));

    //sndbuff has a very non-uniform layout, so we will simply things by moving
    //contiguously through it, and jumping around in in.
    //in = [my_kx->width][my_ky->width][nz]
    //sndbuff = [p][my_kx->width][my_ky->width][pz]
    complex double * pisbuff = sndbuff;
    complex double *  piin = in;

    int offset = 0;
    for(i = 0; i < vsize; i++)
    {
        for(j = 0; j < my_kx->width; j++)
        {
            for(k = 0; k < my_ky->width; k++)
            {
                piin = in + offset + k * nz + j * my_ky->width * nz;
                memcpy(pisbuff, piin, all_z[i].width * sizeof(complex double));
                pisbuff += all_z[i].width;
            }
        }
        offset += all_z[i].width;
    }

    //For one of the processors, the information needed is not contiguous because
    //it is interrupted by wavelengths we wish to discard for dealiasing.
    int dProc = -1;
    int nlow = -1;
    int nhigh = -1;
    for(i = 0; i < vsize; i++)
    {
        if(all_kx[i].min <= dealias_kx.min && all_kx[i].max >= dealias_kx.min)
        {
            dProc = i;
            nlow = dealias_kx.min - all_kx[i].min;
            nhigh = all_kx[i].width - nlow;
        }
    }

    rcnt[0] = 2 * all_kx[0].width * my_ky->width * my_z->width;
    scnt[0] = 2 * my_kx->width * my_ky->width * all_z[0].width;
    sdisp[0] = 0;
    rdisp[0] = 0;
    for(i = 1; i < vsize; i++)
    {
        rcnt[i] = 2 * all_kx[i].width * my_ky->width * my_z->width;
        scnt[i] = 2 * my_kx->width * my_ky->width * all_z[i].width;
        sdisp[i] = sdisp[i-1] + scnt[i-1];
        rdisp[i] = rdisp[i-1] + rcnt[i-1];

        //If we are the proc after dProc, we actually have a bunch of empty
        //wave modes to skip over...
        if(i == dProc+1)
        {
            rdisp[i] += 2 * dealias_kx.width * my_ky->width * my_z->width;
        }
    }

    MPI_Alltoallv(sndbuff, scnt, sdisp, MPI_DOUBLE, out, rcnt, rdisp, MPI_DOUBLE, vcomm);

    //now move the data so that we have the dealiased wavelengths back
    memmove(out + (dealias_kx.max+1) * my_ky->width * my_z->width, out + dealias_kx.min * my_ky->width * my_z->width, nhigh * my_ky->width * my_z->width * sizeof(complex double));
    memset(out + dealias_kx.min * my_ky->width * my_z->width, 0, dealias_kx.width * my_ky->width * my_z->width * sizeof(complex double));

    free(sndbuff);
    free(rcnt);
    free(rdisp);
    free(scnt);
    free(sdisp);
}

void testfft1()
{
    int i,j,k,l;
    double  * start = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    double * finish = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    complex double  * comp = (complex double *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex double));

    int len;
    int * ks;
    if(grank == 0)
    {

        len = rand() % 30+5;
        ks = (int*)malloc(len*3*sizeof(int));

        for(i = 0; i < len; i++)
        {
            ks[3*i] = rand()%dealias_kx.min;
            ks[3*i+1] = rand()%dealias_ky.min;
            ks[3*i+2] = rand()%dealias_kz.min;
        }
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, ccomm);

    if(grank != 0)
        ks = (int*)malloc(len*3*sizeof(int));

    MPI_Bcast(ks, len*3, MPI_INT, 0, ccomm);

    generateFunc(ks, len, (double*)start);

    fft1_forward((double*)start, (complex double*)comp);

    int match;
    int k1,k2,k3;
    int index;
    int mcount = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < ndkz; k++)
            {
                index = k + j*ndkz + i*my_ky->width * ndkz;
                double abs = fabs(creal(comp[index])) + fabs(cimag(comp[index]));
                if(abs > 1e-8)
                {
                    //fprintf(stderr, "found: %d %d %d\n", i + my_kx->min, j + my_ky->min, k);
                    match = 0;
                    k1 = i + my_kx->min;
                    if(k1 >= dealias_kx.min)
                        k1 = ndkx - k1;
                    k2 = j + my_ky->min;
                    k3 = k;
                    if(k3 >= dealias_kz.min)
                        k3 = ndkz - k3;
                    for(l = 0; l < len; l++)
                    {
                        if(ks[3*l] == k1 && ks[3*l+1] == k2 && ks[3*l+2] == k3)
                        {
                            match = 1;
                            mcount++;
                            break;
                        }
                    }
                    if(!match)
                    {
                        fprintf(stderr, "NO MATCH: ");
                        fprintf(stderr,"%d %d %d %d %e + %e i\n", i + my_kx->min, j + my_ky->min, k, grank, creal(comp[index]), cimag(comp[index]));
                    }
                }
            }
        }
    }
    int total;
    MPI_Reduce(&mcount, &total,1, MPI_INT, MPI_SUM, 0, ccomm);
    if(grank == 0)
        fprintf(stderr, "Matched %d out of %d\n", total, len);

    fft1_backward((complex double*)comp, (double*)finish);

    double err;
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_x->width; j++)
        {
            for(k = 0; k < ny; k++)
            {
                int index = k + j * ny + i * my_x->width * ny;
                err = fabs(finish[index] - start[index]);
                if(err > 1e-10)
                    fprintf(stderr, "Problem found! %g %g %g\n", err, finish[index], start[index]);
                //fprintf(stderr, "%d: %d %d %d %g\n", grank, i, j, k, finish[i][j][k]);
            }
        }
    }
}

void testfft2()
{
    int i,j,k,l;
    double  * start = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    double * finish = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    complex double  * comp = (complex double *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex double));

    int len;
    int * ks;
    if(grank == 0)
    {

        len = rand() % 30+5;
        ks = (int*)malloc(len*3*sizeof(int));

        for(i = 0; i < len; i++)
        {
            ks[3*i] = rand()%dealias_kx.min;
            ks[3*i+1] = rand()%dealias_ky.min;
            ks[3*i+2] = rand()%dealias_kz.min;
        }
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, ccomm);

    if(grank != 0)
        ks = (int*)malloc(len*3*sizeof(int));

    MPI_Bcast(ks, len*3, MPI_INT, 0, ccomm);

    generateFunc(ks, len, (double*)start);

    fft2_forward((double*)start, (complex double*)comp);

    int match;
    int k1,k2,k3;
    int index;
    int mcount = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < ndkz; k++)
            {
                index = k + j*ndkz + i*my_ky->width * ndkz;
                double abs = fabs(creal(comp[index])) + fabs(cimag(comp[index]));
                if(abs > 1e-8)
                {
                    match = 0;
                    k1 = i + my_kx->min;
                    if(k1 >= dealias_kx.min)
                        k1 = ndkx - k1;
                    k2 = j + my_ky->min;
                    k3 = k;
                    if(k3 >= dealias_kz.min)
                        k3 = ndkz - k3;
                    for(l = 0; l < len; l++)
                    {
                        if(ks[3*l] == k1 && ks[3*l+1] == k2 && ks[3*l+2] == k3)
                        {
                            match = 1;
                            mcount++;
                            break;
                        }
                    }
                    if(!match)
                    {
                        fprintf(stderr, "NO MATCH: ");
                        fprintf(stderr,"%d %d %d %d %e + %e i\n", i + my_kx->min, j + my_ky->min, k, grank, creal(comp[index]), cimag(comp[index]));
                    }
                }
            }
        }
    }
    int total;
    MPI_Reduce(&mcount, &total,1, MPI_INT, MPI_SUM, 0, ccomm);
    if(grank == 0)
        fprintf(stderr, "Matched %d out of %d\n", total, len);

    fft2_backward((complex double*)comp, (double*)finish);

    double err;
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_x->width; j++)
        {
            for(k = 0; k < ny; k++)
            {
                int index = k + j * ny + i * my_x->width * ny;
                err = fabs(finish[index] - start[index]);
                if(err > 1e-10)
                {
                    fprintf(stderr, "%d Problem found! %g %g %g\n", grank, err, finish[index], start[index]);
                    return;
                }
            }
        }
    }
}

void generateFunc(int* ks, int len, double* out)
{
    int i,j,k,l;
    int * piks;
    double * piout = out;
    double x,y,z;
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_x->width; j++)
        {
            for(k = 0; k < ny; k++)
            {   
                x = 2 * PI * ((double)j + my_x->min) / (double)nx;
                y = 2 * PI * ((double)k) / (double)ny;
                z = 2 * PI * ((double)i + my_z->min) / (double)nz;

                piks = ks;
                (*piout) = 0;
                for(l = 0; l < len; l++)
                {
                    (*piout) += sin(piks[0] * x + piks[1]*y + piks[2]*z + PI / 4.);
                    piks += 3;
                }
                piout++;
            }
        }
    }
}

void repeatfft1(int count)
{
    int i;
    double  * start = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    double * finish = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    complex double  * comp = (complex double *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex double));

    int len;
    int * ks;
    if(grank == 0)
    {

        len = rand() % 30+5;
        ks = (int*)malloc(len*3*sizeof(int));

        for(i = 0; i < len; i++)
        {
            ks[3*i] = rand()%dealias_kx.min;
            ks[3*i+1] = rand()%dealias_ky.min;
            ks[3*i+2] = rand()%dealias_kz.min;
        }
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, ccomm);

    if(grank != 0)
        ks = (int*)malloc(len*3*sizeof(int));

    MPI_Bcast(ks, len*3, MPI_INT, 0, ccomm);

    generateFunc(ks, len, (double*)start);

    for(i = 0; i < count; i++)
    {
        fft1_forward((double*)start, (complex double*)comp);
        fft1_backward((complex double*)comp, (double*)finish);
    }
    fftw_cleanup();
}

void repeatfft2(int count)
{
    int i;
    double  * start = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    double * finish = (double *)malloc(my_z->width * my_x->width * ny * sizeof(double));
    complex double  * comp = (complex double *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex double));

    int len;
    int * ks;
    if(grank == 0)
    {

        len = rand() % 30+5;
        ks = (int*)malloc(len*3*sizeof(int));

        for(i = 0; i < len; i++)
        {
            ks[3*i] = rand()%dealias_kx.min;
            ks[3*i+1] = rand()%dealias_ky.min;
            ks[3*i+2] = rand()%dealias_kz.min;
        }
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, ccomm);

    if(grank != 0)
        ks = (int*)malloc(len*3*sizeof(int));

    MPI_Bcast(ks, len*3, MPI_INT, 0, ccomm);

    generateFunc(ks, len, (double*)start);

    for(i = 0; i < count; i++)
    {
        fft2_forward((double*)start, (complex double*)comp);
        fft2_backward((complex double*)comp, (double*)finish);
    }
    fftw_cleanup();
}

void fftForward(p_field f)
{
    if(whichfft == FFT1)
        fft1_forward(f->spatial, f->spectral);
    else
        fft2_forward(f->spatial, f->spectral);
}

void fftBackward(p_field f)
{
    if(whichfft == FFT1)
        fft1_backward(f->spectral, f->spatial);
    else
        fft2_backward(f->spectral, f->spatial);
}

