/*
 * Copywrite 2013 Benjamin Byington
 *
 * This file is part of the IMHD software package
 * 
 * IMHD is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
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

/********************* 
 * This file controls all the state variables that persist through the
 * simulation and stores the current values of all state variables in the 
 * problem  
 *********************/

#include "Field.h"
#include "Precision.h"

#ifndef _STATE_H
#define	_STATE_H

/*
 * Init and cleanup routines for the state variables.  Should be called once
 * each, before and after calculations.
 */
void initState();
void finalizeState();

//These are our state variables in this system.
extern p_componentVar B;
extern p_componentVar u;
extern p_field T;

//field to be used for boundary hyper diffusion term.  Experimental and not
//fully tested!
extern p_field hyper;
extern p_field hyperWork;

//convenience variable used in status file outputs.
extern PRECISION maxVel[3];

//If we have time independent forcings, they get loaded into here.
extern p_field forceField;
extern p_field magForceField;


#endif	/* _STATE_H */

