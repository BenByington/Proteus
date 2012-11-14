/* 
 * File:   FFTWrapper.h
 * Author: Ben
 *
 * Created on November 18, 2010, 1:24 PM
 */

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

