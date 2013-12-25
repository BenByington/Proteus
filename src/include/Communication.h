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
 * File:   Communication.h
 * Author: Ben
 *
 *
 * 
 * Created on February 18, 2010, 1:15 PM
 */


//TODO: Verify that COMM methods should only be called by compute nodes
#ifndef _COMMUNICATION_H
#define	_COMMUNICATION_H

#include "Field.h"
#include "Precision.h"

#include <mpi.h>
#include <math.h>

#define FFT1 1
#define FFT2 2

void com_init(int measure);
void com_finalize();

void fftForward(p_field);
void fftBackward(p_field);


#endif	/* _COMMUNICATION_H */

