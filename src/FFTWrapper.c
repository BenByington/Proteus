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

#include "FFTWrapper.h"


FFT_PLAN fft_plan_r2c(int rank, const int *n, int howmany,
                          PRECISION *in, const int *inembed,
                          int istride, int idist,
                          FFT_COMPLEX *out, const int *onembed,
                          int ostride, int odist,
                          unsigned flags)
{
    #ifdef FP
    return fftwf_plan_many_dft_r2c(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags);
    #else
    return fftw_plan_many_dft_r2c(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags);
    #endif
}

FFT_PLAN fft_plan_c2c(int rank, const int *n, int howmany,
                      FFT_COMPLEX *in, const int *inembed,
                      int istride, int idist,
                      FFT_COMPLEX *out, const int *onembed,
                      int ostride, int odist,
                      int sign, unsigned flags)
{

    #ifdef FP
    return fftwf_plan_many_dft(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags);
    #else
    return fftw_plan_many_dft(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags);
    #endif

}

FFT_PLAN fft_plan_c2r(int rank, const int *n, int howmany,
                      FFT_COMPLEX *in, const int *inembed,
                      int istride, int idist,
                      PRECISION *out, const int *onembed,
                      int ostride, int odist,
                      unsigned flags)
{
    #ifdef FP
    return fftwf_plan_many_dft_c2r(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags);
    #else
    return fftw_plan_many_dft_c2r(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags);
    #endif
}

void * fft_malloc(size_t size)
{
    #ifdef FP
    return fftwf_malloc(size);
    #else
    return fftw_malloc(size);
    #endif
}

void fft_free(void * data)
{
    #ifdef FP
    fftwf_free(data);
    #else
    fftw_free(data);
    #endif
}

void fft_execute_r2c(FFT_PLAN plan, PRECISION * in, FFT_COMPLEX * out)
{
    #ifdef FP
    fftwf_execute_dft_r2c(plan, in, out);
    #else
    fftw_execute_dft_r2c(plan, in, out);
    #endif
}

void fft_execute_c2c(FFT_PLAN plan, FFT_COMPLEX * in, FFT_COMPLEX * out)
{
    #ifdef FP
    fftwf_execute_dft(plan, in, out);
    #else
    fftw_execute_dft(plan, in, out);
    #endif
}

void fft_execute_c2r(FFT_PLAN plan, FFT_COMPLEX * in, PRECISION * out)
{
    #ifdef FP
    fftwf_execute_dft_c2r(plan, in, out);
    #else
    fftw_execute_dft_c2r(plan, in, out);
    #endif
}

