/* 
 * File:   Precision.h
 * Author: Ben
 *
 * Created on November 17, 2010, 3:34 PM
 */

#ifndef _PRECISION_H
#define	_PRECISION_H



#ifdef FP

#define PRECISION float
#define MPI_PRECISION MPI_FLOAT
#define FFT_PLAN fftwf_plan
#define FFT_COMPLEX fftwf_complex

#else

#define PRECISION double
#define MPI_PRECISION MPI_DOUBLE
#define FFT_PLAN fftw_plan
#define FFT_COMPLEX fftw_complex

#endif


#endif	/* _PRECISION_H */

