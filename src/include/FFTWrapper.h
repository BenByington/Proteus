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

/***************************
 * This signature for FFTW routines depends on the precision the data is 
 * stored in.  This is a lightweight interface designed to hide that dependence,
 * as the precision this code runs in is determined at compile time. 
 ***************************/

#ifndef _FFTWRAPPER_H
#define	_FFTWRAPPER_H

#include "Precision.h"
#include <fftw3.h>

#ifdef	__cplusplus
extern "C" {
#endif

    void * fft_malloc(size_t in);
    void fft_free(void * in);

    void fft_execute_r2c(FFT_PLAN plan, PRECISION * in, FFT_COMPLEX * out);
    void fft_execute_c2c(FFT_PLAN plan, FFT_COMPLEX * in, FFT_COMPLEX * out);
    void fft_execute_c2r(FFT_PLAN plan, FFT_COMPLEX * in, PRECISION * out);

    FFT_PLAN fft_plan_r2c(int rank, const int *n, int howmany,
                          PRECISION *in, const int *inembed,
                          int istride, int idist,
                          FFT_COMPLEX *out, const int *onembed,
                          int ostride, int odist,
                          unsigned flags);

    FFT_PLAN fft_plan_c2c(int rank, const int *n, int howmany,
                          FFT_COMPLEX *in, const int *inembed,
                          int istride, int idist,
                          FFT_COMPLEX *out, const int *onembed,
                          int ostride, int odist,
                          int sign, unsigned flags);

    FFT_PLAN fft_plan_c2r(int rank, const int *n, int howmany,
                          FFT_COMPLEX *in, const int *inembed,
                          int istride, int idist,
                          PRECISION *out, const int *onembed,
                          int ostride, int odist,
                          unsigned flags);


#ifdef	__cplusplus
}
#endif

#endif	/* _FFTWRAPPER_H */

