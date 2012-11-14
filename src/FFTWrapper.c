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