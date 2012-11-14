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
void comm_finalize();

void fftForward(p_field);
void fftBackward(p_field);


#endif	/* _COMMUNICATION_H */

