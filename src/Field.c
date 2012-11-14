/* 
 * File:   Field.cpp
 * Author: Ben
 * 
 * Created on February 12, 2010, 2:06 PM
 */

#include "Field.h"
#include "Environment.h"

#include <fftw3.h>
#include <stdlib.h>

void allocateSpectral(p_field f)
{
    f->spectral = (complex double*)fftw_malloc(my_ky->width * my_kx->width * ndkz * sizeof(complex double));
}

void allocateSpatial(p_field f)
{
    f->spatial = (double*)fftw_malloc(my_x->width * my_z->width * ny * sizeof(double));
}

void allocateForce(p_field f)
{
    f->force1 = (complex double*)fftw_malloc(spectralCount * sizeof(complex double));
    f->force2 = (complex double*)fftw_malloc(spectralCount * sizeof(complex double));
    f->force3 = (complex double*)fftw_malloc(spectralCount * sizeof(complex double));
}

void eraseSpatial(p_field f)
{
    if(f->spatial)
    {
        fftw_free(f->spatial);
        f->spatial = 0;
    }
}

void eraseSpectral(p_field f)
{
    if(f->spectral)
    {
        fftw_free(f->spectral);
        f->spectral = 0;
    }
}

void eraseForce(p_field f)
{
    if(f->force1)
    {
        fftw_free(f->force1);
        f->force1 = 0;
    }

    if(f->force2)
    {
        fftw_free(f->force2);
        f->force2 = 0;
    }

    if(f->force3)
    {
        fftw_free(f->force3);
        f->force3 = 0;
    }
}

p_solenoid newSolenoid()
{
   p_solenoid ret = (p_solenoid)malloc(sizeof(solenoid));

   ret->poloidal = (p_field)malloc(sizeof(field));
   ret->poloidal->spatial = 0;
   ret->poloidal->spectral = 0;
   allocateSpectral(ret->poloidal);
   allocateForce(ret->poloidal);

   ret->toroidal = (p_field)malloc(sizeof(field));
   ret->toroidal->spatial = 0;
   ret->toroidal->spectral = 0;
   allocateSpectral(ret->toroidal);
   allocateForce(ret->toroidal);

   ret->mean_z = 0;

   ret->mean_x = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_y = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_xf1 = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_yf1 = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_xf2 = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_yf2 = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_xf3 = (complex double*)malloc(ndkz * sizeof(complex double));
   ret->mean_yf3 = (complex double*)malloc(ndkz * sizeof(complex double));

   return ret;
}

void deleteSolenoid(p_solenoid * pv)
{
    eraseSpectral((*pv)->poloidal);
    eraseForce((*pv)->poloidal);
    free((*pv)->poloidal);

    eraseSpectral((*pv)->toroidal);
    eraseForce((*pv)->toroidal);
    free((*pv)->toroidal);

    free((*pv)->mean_x);
    free((*pv)->mean_y);
    free((*pv)->mean_xf1);
    free((*pv)->mean_yf1);
    free((*pv)->mean_xf2);
    free((*pv)->mean_yf2);
    free((*pv)->mean_xf3);
    free((*pv)->mean_yf3);
    
    free(*pv);
    *pv = 0;
}

p_vector newVector(int alloc)
{
    p_vector ret = (p_vector)malloc(sizeof(vector));

    ret->x = (p_field)malloc(sizeof(field));
    if(alloc & SPAT)
        allocateSpatial(ret->x);
    else
        ret->x->spatial = 0;
    if(alloc & SPEC)
        allocateSpectral(ret->x);
    else
        ret->x->spectral = 0;

    ret->y = (p_field)malloc(sizeof(field));
    if(alloc & SPAT)
        allocateSpatial(ret->y);
    else
        ret->y->spatial = 0;
    if(alloc & SPEC)
        allocateSpectral(ret->y);
    else
        ret->y->spectral = 0;

    ret->z = (p_field)malloc(sizeof(field));
    if(alloc & SPAT)
        allocateSpatial(ret->z);
    else
        ret->z->spatial = 0;
    if(alloc & SPEC)
        allocateSpectral(ret->z);
    else
        ret->z->spectral = 0;

    return ret;
}

void deleteVector(p_vector * pv)
{
    eraseSpatial((*pv)->x);
    eraseSpectral((*pv)->x);
    free((*pv)->x);

    eraseSpatial((*pv)->y);
    eraseSpectral((*pv)->y);
    free((*pv)->y);

    eraseSpatial((*pv)->z);
    eraseSpectral((*pv)->z);
    free((*pv)->z);

    free(*pv);
    *pv = 0;
}

p_componentVar newComponentVar()
{
    p_componentVar ret = (p_componentVar)malloc(sizeof(componentVar));

    ret->sol = newSolenoid();
    ret->vec = newVector(SPEC | SPAT);

    return ret;
}

void deleteComponentVar(p_componentVar * pc)
{
    deleteVector(&(*pc)->vec);
    deleteSolenoid(&(*pc)->sol);

    (*pc) = 0;
}
