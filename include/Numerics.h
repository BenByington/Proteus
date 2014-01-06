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

#include "Precision.h"
#include "Field.h"
#include "Environment.h"

/*
 * This is a collection of routines that know how to do some of the more
 * common mathematical operators on the data variables as defined in this
 * program.
 */

/*
 * Unit test to ensure that the poloidal and toroidal decomposition 
 * function as expected.  It is currently called at the beginning of
 * execution, but really should only be called when we explicitly have reason
 * to verify things.
 */
void testPT();

/*
 * Derivatives are done in spectral space, where they reduce to a multiplication
 * by wavenumber, and this factor is provided by these routines.
 */
std::complex<PRECISION> dxFactor(int i);
std::complex<PRECISION> dyFactor(int i);
std::complex<PRECISION> dzFactor(int i);

/*
 * These methods deal with the poloidal and toroidal decomposition.  There are
 * two decomposition routines, as sometimes we have an incompressible force and
 * sometimes we have the curl of an incompressible force.  For both decompose
 * routines, the force flag indicates if we should store the result in the
 * force vectors of the appropriate variables.
 */
void decomposeSolenoidal(p_solenoid s, p_vector v, int force);
void decomposeCurlSolenoidal(p_solenoid s, p_vector v, int force);
void recomposeSolenoidal(p_solenoid s, p_vector v);

/*
 * Standard divergence, gradient ad curl.  First argument is the field the 
 * argument is applied to and the second argument is where the result is to be
 * stored.
 */
void curl(p_vector in, p_vector out);
void gradient(p_field in, p_vector out);
void divergence(p_vector in, p_field out);

/*
 * Standard dot and cross product.  The first two arguments are the vectors
 * being operated on, and the result is stored in the third.
 */
void dotProduct(p_vector one, p_vector two, p_field out);
void crossProduct(p_vector one, p_vector two, p_vector out);

/*
 * Routines to take the derivative on a single spectral field.  The first
 * argument is the field we take the derivative of and the second is the
 * resultant.  The arithmetic argument is a flag dictating how the results are
 * stored:
 * 
 * 0 -- out = derivative
 * 1 -- out += derivative
 * 2 -- out -= derivative
 */
void partialX(std::complex<PRECISION> * in, std::complex<PRECISION> * out, int arithmetic);
void partialY(std::complex<PRECISION> * in, std::complex<PRECISION> * out, int arithmetic);
void partialZ(std::complex<PRECISION> * in, std::complex<PRECISION> * out, int arithmetic);

/*
 * More basic operations, again with the first two arguments being the subject
 * of the operation with the third argument holding the result
 */
void multiply(PRECISION * one, PRECISION * two, PRECISION * out);
void plusEq(std::complex<PRECISION> * one, std::complex<PRECISION> * two);
void minusEq(std::complex<PRECISION> * one, std::complex<PRECISION> * two);

/*
 * Experimental and unfinished routines.  These ones work, but they are 
 * mostly useless without the sanitation routines.  
 * 
 * These routines are used to shift a given spectral field f by a given
 * displacement d.
 */
void shiftField(displacement d, std::complex<PRECISION> * f);
void shiftAvg(displacement d, std::complex<PRECISION> * f);

//derivatives for a field
void laplacian(std::complex<PRECISION> * in, std::complex<PRECISION> * out, int add, PRECISION factor);
void hyperDiff(std::complex<PRECISION> * in, std::complex<PRECISION> * out, int add, PRECISION factor);

/*
 * Experimental routine that does not work as intended.  This was originally
 * intended as a way to get an infinite domain, where magnetic field could
 * rise indefinitely.  The state variables would be shifted to keep the
 * magnetic structure centered.  That operation involves "dirty" values being
 * shifted through the periodic boundaries, so this routine was designed to 
 * wipe those out while still maintaining the divergence free constraint.  It
 * fails miserably...
 */
void killBoundaries(PRECISION * in, std::complex<PRECISION> * out, int add, PRECISION factor);

