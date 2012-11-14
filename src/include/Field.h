/* 
 * File:   Field.h
 * Author: Ben
 *
 * Created on February 12, 2010, 2:06 PM
 */

#ifndef _FIELD_H
#define	_FIELD_H

#include "Precision.h"
#include <complex.h>

#define SPEC 1
#define SPAT 2

typedef struct
{
    PRECISION * spatial;
    complex PRECISION * spectral;
    complex PRECISION * force1;
    complex PRECISION * force2;
    complex PRECISION * force3;
}*p_field,field;

typedef struct
{
    p_field poloidal;
    p_field toroidal;
    complex PRECISION * mean_x;
    complex PRECISION * mean_y;
    complex PRECISION * mean_xf1;
    complex PRECISION * mean_yf1;
    complex PRECISION * mean_xf2;
    complex PRECISION * mean_yf2;
    complex PRECISION * mean_xf3;
    complex PRECISION * mean_yf3;
    complex PRECISION mean_z;
}*p_solenoid, solenoid;

typedef struct
{
    p_field x;
    p_field y;
    p_field z;
}*p_vector, vector;

typedef struct
{
    p_vector vec;
    p_solenoid sol;
}*p_componentVar, componentVar;


p_solenoid newSolenoid();
void deleteSolenoid(p_solenoid * ps);

p_vector newVector(int alloc);
void deleteVector(p_vector * pv);

p_componentVar newComponentVar();
void deleteComponentVar(p_componentVar * pc);

void allocateSpectral(field * f);
void allocateSpatial(field * f);
void allocateForce(p_field f);
void eraseSpectral(p_field f);
void eraseSpatial(p_field f);
void eraseForce(p_field f);

#endif	/* _FIELD_H */

