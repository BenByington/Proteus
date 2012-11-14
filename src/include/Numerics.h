/* 
 * File:   Numerics.h
 * Author: Ben
 *
 * Created on April 5, 2010, 10:44 AM
 */

#include "Field.h"

void testPT();

inline complex double dxFactor(int i);
inline complex double dyFactor(int i);
inline complex double dzFactor(int i);

void decomposeSolenoidal(p_solenoid s, p_vector v, int force);
void decomposeCurlSolenoidal(p_solenoid s, p_vector v, int force);
void recomposeSolenoidal(p_solenoid s, p_vector v);

inline void curl(p_vector in, p_vector out);
inline void gradient(p_field in, p_vector out);
inline void divergence(p_vector in, p_field out);
inline void dotProduct(p_vector one, p_vector two, p_field out);
inline void crossProduct(p_vector one, p_vector two, p_vector out);

inline void partialX(complex double * in, complex double * out, int arithmetic);
inline void partialY(complex double * in, complex double * out, int arithmetic);
inline void partialZ(complex double * in, complex double * out, int arithmetic);
inline void multiply(double * one, double * two, double * out);
inline void plusEq(complex double * one, complex double * two);
inline void minusEq(complex double * one, complex double * two);

//derivatives for a field
inline void laplacian(complex double * in, complex double * out, int add, double factor);

