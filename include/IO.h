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

#ifndef _IO_H
#define	_IO_H

#include "Precision.h"
#include <string>
#include "Field.h"

/*
 * Basic test routine to make sure we can write a file and read it back in
 * without corruption.
 */
void testIO();

/*
 * Takes care of some background initialization required for the IO functions.
 * Should only be called once, at the start of execution.
 */
void initIO();

/*
 * Cleans up after the initIO() function.  Should be called once at program
 * termination
 */
void finalizeIO();
 
/*
 * Read-write bulk outputs.  These functions behaves differently depending on if
 * it is called by an IO node or a Compute node. 
 * 
 * For a compute node:
 *   f    :  Field we wish to write to file
 *   name :  <ignored>
 * For an IO node:
 *   f    :  <ignored>
 *   name :  path to the file we wish to write, relative to the directory that
 *           the simulation is running in.
 * 
 * Note: This routine involves MPI collectives, so EVERY processor must call it 
 *       in order to avoid deadlocks. 
 */
void writeSpatial(p_field f, char name[]);
void readSpatial(p_field f, char name[]);

/*
 * Read-write checkpoints.  The location is determined automatically, and only
 * compute nodes should call these routines.
 */
void writeCheckpoint();
void readCheckpoint();

/*
 * Entry point for IO operations.  Call this routine once per iteration, and
 * it will automatically perform the various types of IO operations as needed.
 */
void performOutput();

#endif	/* _IO_H */

