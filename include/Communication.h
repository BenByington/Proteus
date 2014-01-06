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

/***********************
 * This file is the workhorse of the pseudo-spectral algorithm.  We will
 * repeatedly take FFTs of our data, and virtually all of the parallelism of
 * this code is wrapped up in here.
 * 
 * The FFTs should be initialized before use and cleaned up at the end of
 * execution.  Usage is simple: You simply call fftFoward and fftBackward
 * on whichever data you wish transformed.  Forward converts spatial variables
 * to spectral variables, and Backwards does the inverse.  
 ***********************/

#ifndef _COMMUNICATION_H
#define	_COMMUNICATION_H

#include "Field.h"
#include "Precision.h"

#include <mpi.h>
#include <math.h>

#define FFT1 1
#define FFT2 2

/*
 * There should be no calls to an fft routine that is not bracketed by these
 * two calls.
 * 
 * There are two possible FFT routines under the hood.  Passing in a nonzero
 * entry for measure will make the program take some time initially to measure
 * which of these two has better performance on this particular machine.  In
 * principle it doesn't seem to matter terribly much.
 */
void com_init(int measure);
void com_finalize();

/*
 * The forward routine converts spatial data to spectral data, and the backwards
 * routine does the opposite.  The p_field argument contains pointers to both
 * types of data, so the results will be stored into the appropriate component
 * of the input.
 */
void fftForward(p_field);
void fftBackward(p_field);


#endif	/* _COMMUNICATION_H */

