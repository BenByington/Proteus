/*
 * Copyright 2013 Benjamin Byington
 *
 * This file is part of the Proteus software package
 * 
 * Proteus is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Proteus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with Proteus.  If not, see <http://www.gnu.org/licenses/>
 */

/**********************
 * The physics module is where the terms in the system of partial differential
 * equations reside. Each equation currently has its own method for evaluating
 * the force each timestep, and each term in the equations is capable of being
 * enabled or disabled via parameters handed in through the configuration file.
 * 
 * Most of the methods here are not publicly accessible.  From the outside,
 * one only needs to initialize and clean up this module, and iteration is
 * controlled by a single method, very creatively named, iterate().
 * 
 * In addition to pure physics, this module is also in charge of calculating 
 * the next time step, which lies on the border between physics and numerics,
 * and also handles the time integration, which actually should be ported
 * over to Numerics.c. 
 **********************/

#ifndef _PHYSICS_H
#define	_PHYSICS_H

#include "Precision.h"
#include "Environment.h"

void iterate();

/*
 * Standard init and cleanup routines.  Only call each once per execution.
 */
void initPhysics();
void finalizePhysics();

/*
 * functions that can be used to find displacement for field shifting if we wish
 * to keep something of interest centered in the domain. Experimental!!!
 */
displacement displacementByCenter();

#endif	/* _PHYSICS_H */

