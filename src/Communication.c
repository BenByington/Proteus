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
 * This file represents the FFT workhorse that is the basis for this pseudo-
 * spectral algorithm.  It is essentially because we take derivatives in
 * spectral coordinates, multiply nonlinear terms in spatial coordinates, and
 * use FFTs to move between the two.  There are a number of design choices that
 * have been made, which result in the code below.
 * 
 * 1.   The FFTs are performed by the FFTW library.  This has the benefit of 
 *      giving us access to a very portable, robust and optimized FFT library.
 *      However, this may prove a hindrance if this code is ported to 
 *      accelerators that use proprietary FFT routines.  
 * 
 * 2.   There is a fundamental choice between doing a distributed FFT algorithm,
 *      and doing FFTs for individual dimensions local in processor and doing 
 *      global transpose operations to change which dimension is contained in
 *      a single processor.  This code takes the latter route.  When this was
 *      originally written, FFTW did have an MPI version, but it was being
 *      phased out because it was apparently performing poorly.  As of Jan 2014
 *      however it seems to be back, so it may be worth re-evaluating this 
 *      decision when one gets the time.
 * 
 * 3.   We must have at least one dimension distributed to make use of MPI
 *      parallelism, but we can chose if we have a slab (2 dimensions contained
 *      in processor) or pencil (one dimension contained in processor) layout.
 *      The pencil layout requires two transpose steps instead of one, but those
 *      transpose steps are within horizontal or vertical compute layers rather
 *      than using a full all-to-all, and with two distributed dimensions we can
 *      scale to larger problems.  
 * 
 * In addition to those design choices, there are two additional items worth
 * knowing to understand the approach taken below:
 * 
 * 4.   Pseudo-spectral methods require a de-aliasing step.  Anytime you
 *      multiply two wave modes together, say modes m and n, the result will be
 *      in modes m+n and m-n.  This means that if you have N modes, and multiply 
 *      mode N by itself, you have contributions to mode 2N which is not
 *      explicitly represented by our discretization.  Additionally, due to
 *      the aliasing that occurs with high wave number, the energy that should
 *      have been deposited in 2N is deposited in other modes where they belong.
 *      To avoid this problem, every time an FFT is taken, the highest 1/3 of
 *      wave modes are discarded.  Any energy that is deposited into these
 *      de-aliased modes is lost -- just another part of the truncation error
 *      involved in any numerical method -- but now calculating the nonlinear
 *      terms no longer deposits energy into unphysical wave modes.  
 * 
 *      Punch line for the above information: After every FFT we are going to
 *      repack our arrays, discarding the 1/3 highest wave modes, which are 
 *      actually stored in the middle third of the array.
 * 
 * 5.   We have to transpose the data to change which dimension is contiguous
 *      in memory for the next FFTW operation, but we have a choice as to how
 *      that transpose takes place.  Regardless we have to shuffle data between
 *      processors, but we can either transpose the local memory manually, or 
 *      else use a more complicated FFTW interface and let the libraries there
 *      take care of it.  We actually here have two implementations, doing it
 *      both ways.  There is an optional measure routine to see which method
 *      is faster for a given machine, though currently it is disabled.  I don't
 *      believe it makes much difference which method is chosen, though this is
 *      something that needs to be verified explicitly.
 */

#include "Communication.h"
#include "mpi.h"
#include "FFTWrapper.h"

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "Environment.h"
#include "Logs/Log.h"

int whichfft;

FFT_PLAN planf1;
FFT_PLAN planf2;
FFT_PLAN planf3;
FFT_PLAN planb1;
FFT_PLAN planb2;
FFT_PLAN planb3;

void initfft1();
void fft1_forward(PRECISION * in, complex PRECISION * out);
void fft1_tpf1(complex PRECISION * in, complex PRECISION * out);
void fft1_tpf2(complex PRECISION * in, complex PRECISION * out);
void fft1_backward(complex PRECISION * in, PRECISION * out);
void fft1_tpb1(complex PRECISION * in, complex PRECISION * out);
void fft1_tpb2(complex PRECISION * in, complex PRECISION * out);

void fft_tpf3(complex PRECISION * in, complex PRECISION * out);
void fft_tpb3(complex PRECISION * in, complex PRECISION * out);

void initfft2();
void fft2_forward(PRECISION * in, complex PRECISION * out);
void fft2_tpf1(complex PRECISION * in, complex PRECISION * out);
void fft2_tpf2(complex PRECISION * in, complex PRECISION * out);
void fft2_backward(complex PRECISION * in, PRECISION * out);
void fft2_tpb1(complex PRECISION * in, complex PRECISION * out);
void fft2_tpb2(complex PRECISION * in, complex PRECISION * out);

void testfft1();
void repeatfft1(int count);
void testfft2();
void repeatfft2(int count);

void generateFunc(int * ks, int len, PRECISION * out);

/*
 * Both types of FFT implemented here have a separate init routine.  If we are
 * to measure, then try them both and leave the best one active.  If not, then
 * call whichever init routine has been specified.
 */
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
            info("Testing fft1\n");
            testfft1();
        }

        debug("Timing fft1\n");
        start = clock();
        repeatfft1(100);
        stop = clock();
        int diff1 = stop - start;

        fftw_cleanup();
        initfft2();
        if(test)
        {
            info("Testing fft2\n");
            testfft2();
        }

        debug("Timing fft2\n");
        start = clock();
        repeatfft2(100);
        stop = clock();

        int diff2 = stop - start;

        if(test)
        {
            info("fft1 finished in %g\n", (PRECISION)diff1 / CLOCKS_PER_SEC);
            info("fft2 finished in %g\n", (PRECISION)diff2 / CLOCKS_PER_SEC);
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
    debug("Initializing fft1...\n");
    PRECISION * real;
    complex PRECISION * comp1;
    complex PRECISION * comp2;

    real = (PRECISION*)fft_malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    comp1 = (complex PRECISION*)fft_malloc(my_z->width * my_x->width * nky * sizeof(complex PRECISION));
    planf1 = fft_plan_r2c(1, &ny, my_x->width * my_z->width, (PRECISION *)real, 0, 1, ny, comp1, 0, 1, nky, FFTW_MEASURE);
    planb1 = fft_plan_c2r(1, &ny, my_x->width * my_z->width, comp1, 0, 1, nky, (PRECISION*)real, 0, 1, ny, FFTW_MEASURE);
    fft_free(real);
    fft_free(comp1);

    comp1 = (complex PRECISION*)fft_malloc(my_z->width * my_ky->width * nx * sizeof(complex PRECISION));
    comp2 = (complex PRECISION*)fft_malloc(my_z->width * my_ky->width * nkx * sizeof(complex PRECISION));
    planf2 = fft_plan_c2c(1, &nx, my_z->width * my_ky->width, comp1, 0, 1, nx, comp2, 0, 1, nkx, FFTW_FORWARD, FFTW_MEASURE);
    planb2 = fft_plan_c2c(1, &nx, my_z->width * my_ky->width, comp2, 0, 1, nkx, comp1, 0, 1, nx, FFTW_BACKWARD, FFTW_MEASURE);
    fft_free(comp1);
    fft_free(comp2);

    comp1 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nz * sizeof(complex PRECISION));
    comp2 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex PRECISION));
    planf3 = fft_plan_c2c(1, &nz, my_kx->width * my_ky->width, comp1, 0, 1, nz, comp2, 0, 1, nkz, FFTW_FORWARD, FFTW_MEASURE);
    planb3 = fft_plan_c2c(1, &nz, my_kx->width * my_ky->width, comp2, 0, 1, nkz, comp1, 0, 1, nz, FFTW_BACKWARD, FFTW_MEASURE);
    fft_free(comp1);
    fft_free(comp2);

    debug("Initialization done\n");
}

void fft1_forward(PRECISION * in, complex PRECISION* out)
{
    int i;

    trace("Begin fftw1 forward transform\n");
    int mySize1 = my_z->width * my_x->width * nky;
    int mySize2 = my_z->width * my_ky->width * nx;
    int mySize3 = my_kx->width * my_ky->width * nz;

    complex PRECISION * comp = (complex PRECISION*)fft_malloc(mySize1*sizeof(complex PRECISION));
    fft_execute_r2c(planf1, in, comp);

    complex PRECISION * comp2 = (complex PRECISION*)fft_malloc(mySize2*sizeof(complex PRECISION));
    fft1_tpf1(comp, comp2);
    fft_free(comp);

    complex PRECISION * comp3 = (complex PRECISION*)fft_malloc(mySize2*sizeof(complex PRECISION));
    fft_execute_c2c(planf2, comp2, comp3);
    fft_free(comp2);

    complex PRECISION * comp4 = (complex PRECISION*)fft_malloc(mySize3*sizeof(complex PRECISION));
    fft1_tpf2(comp3, comp4);
    fft_free(comp3);

    complex PRECISION * comp5 = (complex PRECISION*)fft_malloc(mySize3*sizeof(complex PRECISION));
    fft_execute_c2c(planf3, comp4, comp5);
    fft_free(comp4);

    fft_tpf3(comp5, out);
    fft_free(comp5);

    int size = my_kx->width * my_ky->width * ndkz;
    PRECISION factor = ny * nx * nz;
    for(i = 0; i < size; i++)
        out[i] /= factor;

    trace("Forward fftw competed\n");
}

void fft1_tpf1(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    int maxSize1 = max_z->width * max_x->width * max_ky->width;
    trace("Starting first transpose for forward fft1\n");
    //in is complex PRECISION[my_z->width][my_x->width][nky]
    //sndbuff and rcvbuff is complex PRECISION[hsize][maxz->width][maxx->width][maxky->width]
    complex PRECISION * sndbuff = (complex PRECISION*)malloc(maxSize1 * hsize * sizeof(complex PRECISION));
    complex PRECISION * rcvbuff = (complex PRECISION*)malloc(maxSize1 * hsize * sizeof(complex PRECISION));

    complex PRECISION * piin = in;
    complex PRECISION * pisbuff = sndbuff;

    trace("Packing arrays for MPI all-to-all\n");
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
                memcpy(pisbuff, piin, all_ky[k].width * sizeof(complex PRECISION));

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

    trace("Sending data over network\n");
    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_PRECISION, rcvbuff, 2 * maxSize1, MPI_PRECISION, hcomm);

    trace("Unpacking data from transfer\n");
    //out is complex PRECISION[myz->width][my_ky->width][nx]
    //rcvbuff is complex PRECISION[hsize][maxz->width][maxx->width][maxky->width]
    //loop overy every element in out and set it to the correct value.  This will
    //perform the transpost
    //TODO look into a more efficient way to do this...
    int xBig;
    int xSmall;
    complex PRECISION * pirbuff = rcvbuff;
    complex PRECISION * piout = out;
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

void fft1_tpf2(complex PRECISION* in, complex PRECISION* out)
{
    trace("Starting second transpose for forward fft1\n");
    int i,j,k;
    int maxSize1 = max_z->width * max_kx->width * max_ky->width;

    //in is complex PRECISION[my_z->width][my_ky->width][nkx]
    //sndbuff and rcvbuff is complex PRECISION[vsize][max_z->width][max_ky->width][max_kx->width]
    complex PRECISION * sndbuff = (complex PRECISION*)malloc(maxSize1 * vsize * sizeof(complex PRECISION));
    complex PRECISION * rcvbuff = (complex PRECISION*)malloc(maxSize1 * vsize * sizeof(complex PRECISION));

    complex PRECISION * piin = in;
    complex PRECISION * pisbuff = sndbuff;

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

    trace("Packing arrays for MPI all-to-all\n");
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
                        memcpy(pisbuff, piin, nlow * sizeof(complex PRECISION));

                    piin += nlow + dealias_kx.width;

                    if(nhigh)
                        memcpy(pisbuff + nlow, piin, nhigh * sizeof(complex PRECISION));
                    piin += nhigh;
                }
                else
                {
                    //our internal pointers start where they need to.  Just do the
                    //memcpy
                    memcpy(pisbuff, piin, all_kx[k].width * sizeof(complex PRECISION));

                    //incriment pointers.  Simple for picomp.  Be careful for
                    //pisbuff
                    piin += all_kx[k].width;
                }
                pisbuff += maxSize1;
            }
        }
    }

    trace("Sending data over network\n");
    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_PRECISION, rcvbuff, 2 * maxSize1, MPI_PRECISION, vcomm);

    trace("Unpacking data from transfer\n");
    //out is complex PRECISION[my_kx->width][my_ky->width][nz]
    //rcvbuff is complex PRECISION[hsize][max_z->width][max_ky->width][max_kx->width]
    //loop overy every element in out and set it to the correct value.  This will
    //perform the transpost
    //TODO look into a more efficient way to do this...
    int zBig;
    int zSmall;
    complex PRECISION * pirbuff = rcvbuff;
    complex PRECISION * piout = out;
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

void fft_tpf3(complex PRECISION * in, complex PRECISION * out)
{
    trace("Performing final dealias for forward transform\n");
    int i,j,k;
    //in is [my_kx][my_ky][nkz]
    //out is [my_kx][my_ky][ndkz]

    complex PRECISION * piin = in;
    complex PRECISION * piout = out;
    //loop over x and y to process contiguous arrays

    int len1 = dealias_kz.min;
    int cut = dealias_kz.width;
    int len2 = nkz - len1 - cut;

    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            memcpy(piout, piin, len1 * sizeof(complex PRECISION));

            piout += len1;
            piin += len1 + cut;

            memcpy(piout, piin, len2*sizeof(complex PRECISION));

            piout += len2;
            piin += len2;
        }
    }
}

void fft1_backward(complex PRECISION* in, PRECISION * out)
{
    trace("Begin fft1 backwards transform\n");
    int mySize1 = my_kx->width * my_ky->width * nz;
    int mySize2 = my_z->width * my_ky->width * nkx;
    int mySize3 = my_z->width * my_x->width * nky;

    complex PRECISION * comp = (complex PRECISION*)fft_malloc(mySize1 * sizeof(complex PRECISION));
    fft_tpb3(in, comp);

    complex PRECISION * comp2 = (complex PRECISION*)fft_malloc(mySize1 * sizeof(complex PRECISION));
    fft_execute_c2c(planb3, comp, comp2);
    fft_free(comp);

    complex PRECISION * comp3 = (complex PRECISION*)fft_malloc(mySize2 * sizeof(complex PRECISION));
    fft1_tpb2(comp2, comp3);
    fft_free(comp2);

    complex PRECISION * comp4 = (complex PRECISION*)fft_malloc(mySize2 * sizeof(complex PRECISION));
    fft_execute_c2c(planb2, comp3, comp4);
    fft_free(comp3);

    complex PRECISION * comp5 = (complex PRECISION*)fft_malloc(mySize3 * sizeof(complex PRECISION));
    fft1_tpb1(comp4, comp5);
    fft_free(comp4);

    fft_execute_c2r(planb1, comp5, out);
    fft_free(comp5);
    
    trace("Inverse FFT completed\n");
}

void fft1_tpb1(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    //in is complex PRECISION[my_z->width][my_ky->width][nx]
    //sndbuff is complex PRECISION[hsize][maxz->width][maxx->width][maxky->width]
    int maxSize1 = max_z->width * max_x->width * max_ky->width;
    complex PRECISION* sndbuff = (complex PRECISION*)malloc(maxSize1 * hsize * sizeof(complex PRECISION));
    complex PRECISION* rcvbuff = (complex PRECISION*)malloc(maxSize1 * hsize * sizeof(complex PRECISION));


    //loop overy every element in in and send it to the correct spot in sndbuff.
    //This will perform the transpose
    //TODO look into a more efficient way to do this...
    int xBig;
    int xSmall;
    complex PRECISION * pisbuff = sndbuff;
    complex PRECISION * piin = in;
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

    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_PRECISION, rcvbuff, 2 * maxSize1, MPI_PRECISION, hcomm);

    complex PRECISION * pirbuff = rcvbuff;
    complex PRECISION * piout = out;
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
                memcpy(piout, pirbuff, all_ky[k].width * sizeof(complex PRECISION));

                piout += all_ky[k].width;
                pirbuff += maxSize1;
            }
            //There are some dealiased wavelengths at the end that we need
            //to put back in as 0's
            memset(piout, 0, dealias_ky.width * sizeof(complex PRECISION));
            piout += dealias_ky.width;
        }
    }

    free(sndbuff);
    free(rcvbuff);
}

void fft1_tpb2(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    int maxSize1 = max_z->width * max_kx->width * max_ky->width;
    complex PRECISION * sndbuff = (complex PRECISION*)malloc(maxSize1 * vsize * sizeof(complex PRECISION));
    complex PRECISION * rcvbuff = (complex PRECISION*)malloc(maxSize1 * vsize * sizeof(complex PRECISION));


    //in is complex PRECISION[my_kx->width][my_ky->width][nz]
    //sndbuff is complex PRECISION[vsize][max_z->width][max_ky->width][max_kx->width]
    //loop overy every element in in and set it to the correct value.  This will
    //perform the transpost
    //TODO look into a more efficient way to do this...
    int zBig;
    int zSmall;
    complex PRECISION * pisbuff = sndbuff;
    complex PRECISION * piin = in;
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

    MPI_Alltoall(sndbuff, 2 * maxSize1, MPI_PRECISION, rcvbuff, 2 * maxSize1, MPI_PRECISION, vcomm);

    complex PRECISION * piout = out;
    complex PRECISION * pirbuff = sndbuff;

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
                        memcpy(piout, pirbuff, nlow * sizeof(complex PRECISION));

                    piout += nlow;
                    memset(piout, 0, dealias_kx.width*sizeof(complex PRECISION));
                    piout += dealias_kx.width;

                    if(nhigh)
                        memcpy(piout, pirbuff + nlow, nhigh * sizeof(complex PRECISION));
                    piout += nhigh;
                }
                else
                {
                    //our internal pointers start where they need to.  Just do the
                    //memcpy
                    memcpy(piout, pirbuff, all_kx[k].width * sizeof(complex PRECISION));

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

void fft_tpb3(complex PRECISION * in, complex PRECISION * out)
{
    int i,j;
    //in is [my_kx][my_ky][ndkz]
    //out is [my_kx][my_ky][nkz]

    complex PRECISION * piin = in;
    complex PRECISION * piout = out;
    //loop over x and y to process contiguous arrays

    int len1 = dealias_kz.min;
    int cut = dealias_kz.width;
    int len2 = nkz - len1 - cut;

    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            memcpy(piout, piin, len1 * sizeof(complex PRECISION));
            piout += len1;
            
            memset(piout, 0, cut * sizeof(complex PRECISION));
            piout += cut;

            piin += len1;
            memcpy(piout, piin, len2*sizeof(complex PRECISION));

            piin += len2;
            piout += len2;
        }
    }
}

/*
 * This is the alternate FFT formulation.  Here we are still responsible for
 * shuffling data between processors and repacking the de-aliased arrays, but
 * we don't have to manually transpose anything, and things can be thought of
 * as contiguous memory copies.
 */
void initfft2()
{
    PRECISION * real;
    complex PRECISION * comp1;
    complex PRECISION * comp2;

    real = (PRECISION*)fft_malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    comp1 = (complex PRECISION*)fft_malloc(nky * my_z->width * my_x->width * sizeof(complex PRECISION));
    planf1 = fft_plan_r2c(1, &ny, my_x->width * my_z->width, (PRECISION *)real, 0, 1, ny, comp1, 0, my_x->width * my_z->width, 1, FFTW_MEASURE);
    planb1 = fft_plan_c2r(1, &ny, my_x->width * my_z->width, comp1, 0, my_x->width * my_z->width, 1, (PRECISION*)real, 0, 1, ny, FFTW_MEASURE);
    fft_free(real);
    fft_free(comp1);

    comp1 = (complex PRECISION*)fft_malloc(my_z->width * my_ky->width * nx * sizeof(complex PRECISION));
    comp2 = (complex PRECISION*)fft_malloc(my_z->width * my_ky->width * nkx * sizeof(complex PRECISION));
    planf2 = fft_plan_c2c(1, &nx, my_z->width * my_ky->width, comp1, 0, 1, nx, comp2, 0, my_z->width * my_ky->width, 1, FFTW_FORWARD, FFTW_MEASURE);
    planb2 = fft_plan_c2c(1, &nx, my_z->width * my_ky->width, comp2, 0, my_z->width * my_ky->width, 1, comp1, 0, 1, nx, FFTW_BACKWARD, FFTW_MEASURE);
    fft_free(comp1);
    fft_free(comp2);

    comp1 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nz * sizeof(complex PRECISION));
    comp2 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex PRECISION));
    planf3 = fft_plan_c2c(1, &nz, my_kx->width * my_ky->width, comp1, 0, 1, nz, comp2, 0, 1, nkz, FFTW_FORWARD, FFTW_MEASURE);
    planb3 = fft_plan_c2c(1, &nz, my_kx->width * my_ky->width, comp2, 0, 1, nkz, comp1, 0, 1, nz, FFTW_BACKWARD, FFTW_MEASURE);
    fft_free(comp1);
    fft_free(comp2);
}

void fft2_forward(PRECISION* in, complex PRECISION* out)
{
    complex PRECISION * comp1 = (complex PRECISION*)fft_malloc(nky * my_z->width * my_x->width * sizeof(complex PRECISION));
    fft_execute_r2c(planf1, in, comp1);

    complex PRECISION * comp2 = (complex PRECISION*)fft_malloc(my_ky->width * my_z->width * nx * sizeof(complex PRECISION));
    fft2_tpf1(comp1, comp2);
    fft_free(comp1);

    complex PRECISION * comp3 = (complex PRECISION*)fft_malloc(nkx * my_ky->width * my_z->width * sizeof(complex PRECISION));
    fft_execute_c2c(planf2, comp2, comp3);
    fft_free(comp2);

    complex PRECISION* comp4 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nz * sizeof(complex PRECISION));
    fft2_tpf2(comp3, comp4);
    fft_free(comp3);

    complex PRECISION * comp5 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex PRECISION));
    fft_execute_c2c(planf3, comp4, comp5);
    fft_free(comp4);

    fft_tpf3(comp5, out);
    fft_free(comp5);

    int i;
    PRECISION factor = ny * nx * nz;
    int size = my_kx->width * my_ky->width * ndkz;
    for(i = 0; i < size; i++)
        out[i] /= factor;
}

void fft2_tpf1(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    complex PRECISION * rcvbuff = (complex PRECISION*)malloc(my_ky->width * my_z->width * nx * sizeof(complex PRECISION));
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

    MPI_Alltoallv(in, scnt, sdisp, MPI_PRECISION, rcvbuff, rcnt, rdisp, MPI_PRECISION, hcomm);

    //rcvbuff has a very non-uniform layout, so we will simply things by moving
    //contiguously through it, and jumping around in out.
    //out = [my_ky->width][my_z->width][nx]
    //rcvbuff = [p][my_ky->width][my_z->width][px]
    complex PRECISION * pirbuff = rcvbuff;
    complex PRECISION *  piout = out;

    int offset = 0;
    for(i = 0; i < hsize; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < my_z->width; k++)
            {
                piout = out + offset + k * nx + j * my_z->width * nx;
                memcpy(piout, pirbuff, all_x[i].width * sizeof(complex PRECISION));
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

void fft2_tpf2(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    complex PRECISION * rcvbuff = (complex PRECISION*)malloc(my_kx->width * my_ky->width * nz * sizeof(complex PRECISION));
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
    memmove(in + dealias_kx.min * my_ky->width * my_z->width, in + (dealias_kx.max+1) * my_ky->width * my_z->width, nhigh * my_ky->width * my_z->width * sizeof(complex PRECISION));

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

    MPI_Alltoallv(in, scnt, sdisp, MPI_PRECISION, rcvbuff, rcnt, rdisp, MPI_PRECISION, vcomm);

    //rcvbuff has a very non-uniform layout, so we will simplify things by moving
    //contiguously through it, and jumping around in out.
    //out = [my_kx->width][my_ky->width][nz]
    //rcvbuff = [p][my_kx->width][my_ky->width][pz]
    complex PRECISION * pirbuff = rcvbuff;
    complex PRECISION *  piout = out;

    int offset = 0;
    for(i = 0; i < vsize; i++)
    {
        for(j = 0; j < my_kx->width; j++)
        {
            for(k = 0; k < my_ky->width; k++)
            {
                piout = out + offset + k * nz + j * my_ky->width * nz;
                memcpy(piout, pirbuff, all_z[i].width * sizeof(complex PRECISION));
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

void fft2_backward(complex PRECISION* in, PRECISION* out)
{
    complex PRECISION * comp = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nkz * sizeof(complex PRECISION));
    fft_tpb3(in, comp);

    complex PRECISION * comp1 = (complex PRECISION*)fft_malloc(my_kx->width * my_ky->width * nkz*sizeof(complex PRECISION));
    fft_execute_c2c(planb3, comp, comp1);
    fft_free(comp);

    complex PRECISION * comp2 = (complex PRECISION*)fft_malloc(nkx * my_ky->width * my_z->width*sizeof(complex PRECISION));
    fft2_tpb2(comp1, comp2);
    fft_free(comp1);

    complex PRECISION * comp3 = (complex PRECISION*)fft_malloc(my_ky->width * my_z->width * nx*sizeof(complex PRECISION));
    fft_execute_c2c(planb2, comp2, comp3);
    fft_free(comp2);
    
    complex PRECISION * comp4 = (complex PRECISION*)fft_malloc(nky * my_z->width * my_x->width*sizeof(complex PRECISION));
    fft2_tpb1(comp3, comp4);
    fft_free(comp3);
    
    fft_execute_c2r(planb1, comp4, out);
    fft_free(comp4);
}

void fft2_tpb1(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    complex PRECISION * sndbuff = (complex PRECISION*)malloc(my_ky->width * my_z->width * nx * sizeof(complex PRECISION));
    //set up the send/receive data structures
    int * scnt = (int*)malloc(hsize * sizeof(int));
    int * sdisp = (int*)malloc(hsize * sizeof(int));
    int * rcnt = (int*)malloc(hsize * sizeof(int));
    int * rdisp = (int*)malloc(hsize * sizeof(int));

    //sndbuff has a very non-uniform layout, so we will simply things by moving
    //contiguously through it, and jumping around in out.
    //out = [my_ky->width][my_z->width][nx]
    //rcvbuff = [p][my_ky->width][my_z->width][px]
    complex PRECISION * pisbuff = sndbuff;
    complex PRECISION *  piin = in;

    int offset = 0;
    for(i = 0; i < hsize; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < my_z->width; k++)
            {
                piin = in + offset + k * nx + j * my_z->width * nx;
                memcpy(pisbuff, piin, all_x[i].width * sizeof(complex PRECISION));
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

    MPI_Alltoallv(sndbuff, scnt, sdisp, MPI_PRECISION, out, rcnt, rdisp, MPI_PRECISION, hcomm);

    //make sure the dealiased wavelengths are 0
    memset(out + dealias_ky.min * my_x->width * my_z->width, 0, dealias_ky.width * my_x->width * my_z->width * sizeof(complex PRECISION));

    free(sndbuff);
    free(rcnt);
    free(rdisp);
    free(scnt);
    free(sdisp);
}

void fft2_tpb2(complex PRECISION* in, complex PRECISION* out)
{
    int i,j,k;
    complex PRECISION * sndbuff = (complex PRECISION*)malloc(my_kx->width * my_ky->width * nz * sizeof(complex PRECISION));
    //set up the send/receive data structures
    int * scnt = (int*)malloc(vsize * sizeof(int));
    int * sdisp = (int*)malloc(vsize * sizeof(int));
    int * rcnt = (int*)malloc(vsize * sizeof(int));
    int * rdisp = (int*)malloc(vsize * sizeof(int));

    //sndbuff has a very non-uniform layout, so we will simply things by moving
    //contiguously through it, and jumping around in in.
    //in = [my_kx->width][my_ky->width][nz]
    //sndbuff = [p][my_kx->width][my_ky->width][pz]
    complex PRECISION * pisbuff = sndbuff;
    complex PRECISION *  piin = in;

    int offset = 0;
    for(i = 0; i < vsize; i++)
    {
        for(j = 0; j < my_kx->width; j++)
        {
            for(k = 0; k < my_ky->width; k++)
            {
                piin = in + offset + k * nz + j * my_ky->width * nz;
                memcpy(pisbuff, piin, all_z[i].width * sizeof(complex PRECISION));
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

    MPI_Alltoallv(sndbuff, scnt, sdisp, MPI_PRECISION, out, rcnt, rdisp, MPI_PRECISION, vcomm);

    //now move the data so that we have the dealiased wavelengths back
    memmove(out + (dealias_kx.max+1) * my_ky->width * my_z->width, out + dealias_kx.min * my_ky->width * my_z->width, nhigh * my_ky->width * my_z->width * sizeof(complex PRECISION));
    memset(out + dealias_kx.min * my_ky->width * my_z->width, 0, dealias_kx.width * my_ky->width * my_z->width * sizeof(complex PRECISION));

    free(sndbuff);
    free(rcnt);
    free(rdisp);
    free(scnt);
    free(sdisp);
}

void testfft1()
{
    int i,j,k,l;
    PRECISION  * start = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    PRECISION * finish = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    complex PRECISION  * comp = (complex PRECISION *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex PRECISION));

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

    generateFunc(ks, len, (PRECISION*)start);

    fft1_forward((PRECISION*)start, (complex PRECISION*)comp);

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
                PRECISION abs = fabs(creal(comp[index])) + fabs(cimag(comp[index]));
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

    fft1_backward((complex PRECISION*)comp, (PRECISION*)finish);

    PRECISION err;
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
    PRECISION  * start = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    PRECISION * finish = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    complex PRECISION  * comp = (complex PRECISION *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex PRECISION));

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

    generateFunc(ks, len, (PRECISION*)start);

    fft2_forward((PRECISION*)start, (complex PRECISION*)comp);

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
                PRECISION abs = fabs(creal(comp[index])) + fabs(cimag(comp[index]));
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

    fft2_backward((complex PRECISION*)comp, (PRECISION*)finish);

    PRECISION err;
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

void generateFunc(int* ks, int len, PRECISION* out)
{
    int i,j,k,l;
    int * piks;
    PRECISION * piout = out;
    PRECISION x,y,z;
    for(i = 0; i < my_z->width; i++)
    {
        for(j = 0; j < my_x->width; j++)
        {
            for(k = 0; k < ny; k++)
            {   
                x = 2 * PI * ((PRECISION)j + my_x->min) / (PRECISION)nx;
                y = 2 * PI * ((PRECISION)k) / (PRECISION)ny;
                z = 2 * PI * ((PRECISION)i + my_z->min) / (PRECISION)nz;

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
    PRECISION  * start = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    PRECISION * finish = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    complex PRECISION  * comp = (complex PRECISION *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex PRECISION));

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

    generateFunc(ks, len, (PRECISION*)start);

    for(i = 0; i < count; i++)
    {
        fft1_forward((PRECISION*)start, (complex PRECISION*)comp);
        fft1_backward((complex PRECISION*)comp, (PRECISION*)finish);
    }
    fftw_cleanup();
}

void repeatfft2(int count)
{
    int i;
    PRECISION  * start = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    PRECISION * finish = (PRECISION *)malloc(my_z->width * my_x->width * ny * sizeof(PRECISION));
    complex PRECISION  * comp = (complex PRECISION *)malloc(my_kx->width * my_ky->width * ndkz *  sizeof(complex PRECISION));

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

    generateFunc(ks, len, (PRECISION*)start);

    for(i = 0; i < count; i++)
    {
        fft2_forward((PRECISION*)start, (complex PRECISION*)comp);
        fft2_backward((complex PRECISION*)comp, (PRECISION*)finish);
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

void com_finalize()
{
    fftw_cleanup();
}

