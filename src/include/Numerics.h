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

