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
 * File:   Numerics.h
 * Author: Ben
 *
 * Created on April 5, 2010, 10:44 AM
 */

#include "Precision.h"
#include "Field.h"
#include "Environment.h"

void testPT();

inline complex PRECISION dxFactor(int i);
inline complex PRECISION dyFactor(int i);
inline complex PRECISION dzFactor(int i);

void decomposeSolenoidal(p_solenoid s, p_vector v, int force);
void decomposeCurlSolenoidal(p_solenoid s, p_vector v, int force);
void recomposeSolenoidal(p_solenoid s, p_vector v);

inline void curl(p_vector in, p_vector out);
inline void gradient(p_field in, p_vector out);
inline void divergence(p_vector in, p_field out);
inline void dotProduct(p_vector one, p_vector two, p_field out);
inline void crossProduct(p_vector one, p_vector two, p_vector out);

inline void partialX(complex PRECISION * in, complex PRECISION * out, int arithmetic);
inline void partialY(complex PRECISION * in, complex PRECISION * out, int arithmetic);
inline void partialZ(complex PRECISION * in, complex PRECISION * out, int arithmetic);
inline void multiply(PRECISION * one, PRECISION * two, PRECISION * out);
inline void plusEq(complex PRECISION * one, complex PRECISION * two);
inline void minusEq(complex PRECISION * one, complex PRECISION * two);

inline void shiftField(displacement d, complex PRECISION * f);
inline void shiftAvg(displacement d, complex PRECISION * f);

//derivatives for a field
inline void laplacian(complex PRECISION * in, complex PRECISION * out, int add, PRECISION factor);
inline void hyperDiff(complex PRECISION * in, complex PRECISION * out, int add, PRECISION factor);
inline void killBoundaries(PRECISION * in, complex PRECISION * out, int add, PRECISION factor);

