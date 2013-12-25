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

