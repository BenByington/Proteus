/*
 * Copywrite 2013 Benjamin Byington
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

#ifndef _LABORDIVISION_H
#define	_LABORDIVISION_H

#include "Precision.h"

/*
 * These methods are used to initialize the code.  They are *only* to be called
 * once, and only from the init method inside Environment.c.  The first method
 * initializes the physical problem, the second initializes the communication
 * groups to be used by MPI, and the last one initializes the domain 
 * decompositions to be used.
 */
void lab_initGeometry();
void lab_initGroups();
void lab_initDistributions();

/*
 * Cleanup routine to be called from main when the program is ready to terminate
 */
void lab_finalize();


#endif	/* _LABORDIVISION_H */

