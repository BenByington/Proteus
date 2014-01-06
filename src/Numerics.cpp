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

#include "Numerics.h"
#include "Logs/Log.h"
#include "Environment.h"
#include "Communication.h"
#include "State.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

/*
 * This is a basic test of our poloida/toroidal decomposition.  We begin with
 * a simple vector field with a known form, pass it through our routines, and
 * make sure that the recompose method really is the inverse of the decompose
 * method.
 */
void testPT()
{
    int i,j,k,l;
    int index = 0;
    PRECISION ampM = 3;
    complex<PRECISION> dkx,dky,dkz;
    p_solenoid s = newSolenoid();
    p_vector v = newVector(SPEC);
    p_vector v2 = newVector(SPEC);

    //set the z field with a bunch of random spectral amplitudes
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < ndkz; k++)
            {
                //bz cannot contain a function solely of z
                if(i + j + grank == 0)
                    v->z->spectral[index] = 0;
                else
                    v->z->spectral[index] = ampM * complex<PRECISION>(((PRECISION)rand()) / RAND_MAX, ((PRECISION)rand()) / RAND_MAX);
                index++;
            }
        }
    }
    
    //set the y direction in similar fashion, though there are constraints now
    index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                
                //functions only of y are not allowed
                if(k + i + vrank == 0)        
                {
                    v->y->spectral[index] = 0;
                }
                //functions of only y and z must cancel out with the corresponding portion
                //from the z field
                else if(i + vrank == 0)
                {
                    v->y->spectral[index] = -dkz * v->z->spectral[index] / dky;
                }
                //we are free to chose everything else.
                else
                {
                    v->y->spectral[index] = ampM * complex<PRECISION>(((PRECISION)rand()) / RAND_MAX, ((PRECISION)rand()) / RAND_MAX);
                }
                index++;
            }
        }
    }
    
    //set the x field
    index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                //no x only field
                if(j + hrank + k == 0)
                {
                    v->x->spectral[index] = 0;
                }
                //we can chose any field of y and z
                else if(i + vrank == 0)
                {
                    v->x->spectral[index] = ampM * complex<PRECISION>(((PRECISION)rand()) / RAND_MAX, ((PRECISION)rand()) / RAND_MAX);
                }
                //everything else must cancel with other fields
                else
                {
                    v->x->spectral[index] = -(dkz * v->z->spectral[index] + dky * v->y->spectral[index])/dkx;
                }
                index++;
            }
        }
    }

    //for(i = 0; i < my_kx->width * my_ky->width *ndkz; i++)
      //  fprintf(stderr,"%g %g\n", __real__ v->x->spectral[i], __imag__ v->x->spectral[i]);


    //Decompose and recompose the vectors.  This should be a unitary transform,
    //with v2 = v
    decomposeSolenoidal(s, v, 0);
    recomposeSolenoidal(s, v2);

    //Loop through the arrays and ensure v = v2 to tolerable precision.
    index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < ndkz; k++)
            {
                PRECISION diff = abs(v2->z->spectral[index] - v->z->spectral[index]);
                if(diff > 1e-10)
                    fprintf(stderr,"diffz: %d %d %d %d %g\n", grank, i + my_kx->min, j + my_ky->min, k, diff);

                diff = abs(v2->y->spectral[index] - v->y->spectral[index]);
                if(diff > 1e-10)
                    fprintf(stderr,"diffy: %d %d %d %d %g\n", grank, i + my_kx->min, j + my_ky->min, k, diff);

                diff = abs(v2->x->spectral[index] - v->x->spectral[index]);
                if(diff > 1e-10)
                    fprintf(stderr,"diffx: %d %d %d %d %g\n", grank, i + my_kx->min, j + my_ky->min, k, diff);
                index++;
            }
        }
    }

}

/*
 * We can write a given divergence free vector A as:
 * 
 * A = curl(T hat z) + curl(curl(P hat z)) + < Ax_m(z), Ay_m(z) , Az_m >
 * 
 * In component form this gives us:
 * 
 * curl(T hat z)       = < dy T, -dx T, 0 >
 * curl(curl(P hat z)) = < dx dz P, dy dz P, -(dx dx + dy dy)P >
 * A                   = < dx dz P + dy T, dy dz P - dx T, -(dx dx + dy dy)P >
 * 
 * This is a linear operator, so all of our wavemodes remain independent. It is
 * easiest to separate the modes by their spatial dependencies, e.g.
 * P = F1(x) + F2(y) + F3(z) + F4(x,y) + F5(x,z) + F6(y,z) + F7(x,y,z)
 * 
 * First we apply the inverse horizontal laplacian to the z component of our 
 * vector A.  This uniquely determines all of P save for F3, which we chose to 
 * set to zero.  
 * 
 * Finding the toroidal field happens in two stages.  Any dependence of T on
 * *just* x will not appear in Ax, and similarly any dependence of T on *just*
 * y will not appear in Ay.  So first we use Ax to solve for F2, F3 and F6,
 * since they cannot be gleaned from Ay.  Ay is then used to find the remaining
 * unknowns.
 * 
 * Note that it is possible, for example, to find F7 from either Ax or Ay.  If
 * a properly divergence free vector was passed in, then either way will result
 * in the same answer.  If a non-divergence free vector is passed in, then this
 * decomposition will not be reversible. In practice this means that the small
 * truncation errors inherent in any discretized code are projected onto a 
 * nearby divergence free vector, and any large divergent components will result
 * in unpredictable behavior. 
 */
void decomposeSolenoidal(p_solenoid s, p_vector v, int force)
{
    //TODO: do some code duplication that will allow us to avoid checking
    //for 0 wave modes every iteration of the loops
    int i,j,k;
    debug("Beginning decomposition of solenoidal field\n");

    if(crank == 0)
    {
        s->mean_z = v->z->spectral[0];
    }

    debug("Calculating Poloidal field\n");

    complex<PRECISION> dkx,dky,dkz;
    complex<PRECISION> * pz = v->z->spectral;
    complex<PRECISION> * py = v->y->spectral;
    complex<PRECISION> * px = v->x->spectral;
    complex<PRECISION> * ppol;
    complex<PRECISION> * ptor;
    complex<PRECISION> * xmean;
    complex<PRECISION> * ymean;

    //Are we storing the result in the forcing fields or in the spectral fields.
    if(force)
    {
        ppol = s->poloidal->force1;
        ptor = s->toroidal->force1;
        xmean = s->mean_xf1;
        ymean = s->mean_yf1;
    }
    else
    {
        ppol = s->poloidal->spectral;
        ptor = s->toroidal->spectral;
        xmean = s->mean_x;
        ymean = s->mean_y;
    }


    int index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                if((i+j+crank)==0)
                    ppol[index] = 0;
                else
                    ppol[index] = -pz[index] / (dkx*dkx + dky*dky);

                index++;
            }
        }
    }

    debug("Calculating Toroidal Field\n");
    index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                if(i+j+crank == 0)
                {
                    ptor[index] = 0;
                    xmean[k] = px[index];
                    ymean[k] = py[index];

                }
                //functions of only y and z come from the x field
                else if(i+vrank==0)
                {
                    ptor[index] = (px[index] - dkx * dkz * ppol[index])/dky;
                }
                //everything else we can find from y field
                else
                {
                    ptor[index] = -(py[index] - dky * dkz * ppol[index])/dkx;
                }
                index++;
            }
        }
    }

}

/*
 * Very similar to the previous function, save that now we are dealing with the
 * curl of a divergence free vector and thus have:
 * 
 * curl(T hat z)       = < dy T, -dx T, 0 >
 * curl(curl(P hat z)) = < dx dz P, dy dz P, -(dx dx + dy dy)P >
 * A                   = < dx dz P + dy T, dy dz P - dx T, -(dx dx + dy dy)P >
 * curl(A)             = < dx dz T - dy del^2 P, dy dz T + dx del^2 P, - (dx dx + dy dy) T> 
 * 
 * We proceed as before, save that the inverse derivatives are now slightly
 * more complicated, and the wavenumbers that were previously used to find P are
 * now used to find T, and vise versa. 
 */
void decomposeCurlSolenoidal(p_solenoid s, p_vector v, int force)
{
    //TODO: do some code duplication that will allow us to avoid checking
    //for 0 wave modes every iteration of the loops
    int i,j,k;
    debug("Beginning decomposition of curled solenoidal field\n");

    if(crank == 0)
    {
        s->mean_z = 0;
    }

    debug("Calculating Toroidal field\n");

    complex<PRECISION> dkx,dky,dkz;
    complex<PRECISION> * pz = v->z->spectral;
    complex<PRECISION> * py = v->y->spectral;
    complex<PRECISION> * px = v->x->spectral;
    complex<PRECISION> * ppol;
    complex<PRECISION> * ptor;
    complex<PRECISION> * xmean;
    complex<PRECISION> * ymean;

    if(force)
    {
        ppol = s->poloidal->force1;
        ptor = s->toroidal->force1;
        xmean = s->mean_xf1;
        ymean = s->mean_yf1;
    }
    else
    {
        ppol = s->poloidal->spectral;
        ptor = s->toroidal->spectral;
        xmean = s->mean_x;
        ymean = s->mean_y;
    }


    int index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                if((i+j+crank)==0)
                    ptor[index] = 0;
                else
                    ptor[index] = -pz[index] / (dkx*dkx + dky*dky);

                index++;
            }
        }
    }

    debug("Calculating Poloidal Field\n");
    index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                if(i+j+crank == 0)
                {
                    ppol[index] = 0;
                    if(k==0)
                    {
                        //something curled wont have a box mean
                        //Note: For most things this wont make a difference, but
                        //this implicitly means that you can't chose a forcing 
                        //that gives a variable a net acceleration upwards, 
                        //though it is unclear why that would be wanted in a
                        //periodic domain.
                        xmean[k] = 0;
                        ymean[k] = 0;
                        s->mean_z = 0;

                    }
                    else
                    {
                        xmean[k] = py[index] / dkz;
                        ymean[k] = -px[index] / dkz;
                    }

                }
                //functions of only y and z come from the x field
                else if(i+vrank==0)
                {
                    ppol[index] = -px[index]/dky/(dkx*dkx + dky*dky + dkz*dkz);
                }
                //everything else we can find from y field
                else
                {
                    ppol[index] = (py[index] - dky * dkz * ptor[index])/dkx/(dkx*dkx + dky*dky + dkz*dkz);
                }
                index++;
            }
        }
    }
}

/*
 * Recomposing is much simpler than decomposing.  Here we apply derivatives and
 * directly use:
 * 
 * A  = < dx dz P + dy T, dy dz P - dx T, -(dx dx + dy dy)P >
 */
void recomposeSolenoidal(p_solenoid s, p_vector v)
{
    int i,j,k;
    complex<PRECISION> dkx,dky,dkz;
    complex<PRECISION> * pz = v->z->spectral;
    complex<PRECISION> * py = v->y->spectral;
    complex<PRECISION> * px = v->x->spectral;
    complex<PRECISION> * ppol = s->poloidal->spectral;
    complex<PRECISION> * ptor = s->toroidal->spectral;
    complex<PRECISION> * xmean = s->mean_x;
    complex<PRECISION> * ymean = s->mean_y;

    int index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);

                if(i + j+ crank == 0)
                {
                    px[index] = xmean[k];
                    py[index] = ymean[k];
                }
                else
                {
                    px[index] = dky * ptor[index] + dkx * dkz * ppol[index];
                    py[index] = -dkx * ptor[index] + dky * dkz * ppol[index];
                }

                pz[index] = -(dkx*dkx + dky*dky)*ppol[index];
                
                index++;
            }
        }
    }

    if(crank == 0)
    {
        v->z->spectral[0] = s->mean_z;
    }
}

/*
 * Here we just apply the laplacian, which is a linear operator so each wavemode
 * is simply multiplied by (dx dx + dy dy + dz dz)*factor.  Since each wavemode
 * is independent, this operation can safely be done in place and both out and
 * in can point to the same memory as long as in is no longer needed. 
 */
void laplacian(complex<PRECISION> * in, complex<PRECISION> * out, int add, PRECISION factor)
{
    trace("Starting Laplacian\n");
    
    int i,j,k;
    int index = 0;
    complex<PRECISION> dkx,dky,dkz;

    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                if(add)
                    out[index] += factor * (dkx * dkx + dky * dky + dkz * dkz)*in[index];
                else
                    out[index] = factor * (dkx * dkx + dky * dky + dkz * dkz)*in[index];

                index++;
            }
        }
    }

    trace("Finished Laplacian\n");
}

/*
 * Experimental hyper diffusion used to try and remove field near the boundaries
 * of the computation without breaking the divergence constraints on solenoidal
 * variables.  The hyperdiffusion coefficient is 0 interior to the domain and
 * smoothely ramps up to 'factor' (tanh) at the boundaries.  
 * 
 * DOES NOT WORK AS DESIRED
 */
void hyperDiff(complex<PRECISION> * in, complex<PRECISION> * out, int add, PRECISION factor)
{
    trace("Starting Hyper Diffusion\n");
    
    int i,j,k;
    int index = 0;
    complex<PRECISION> dkx,dky,dkz;
    complex<PRECISION> deriv;
    
    complex<PRECISION> * spect = hyperWork->spectral;
    PRECISION * spat = hyperWork->spatial;

    //Calculate hyper diffusion as if with a constant coefficient
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                deriv = (dkx * dkx + dky * dky + dkz * dkz);
                spect[index] = factor * deriv *in[index];
                index++;
            }
        }
    }
    
    //Transform things to spatial and back, modulating things so that our
    //force vanishes near the interior of the domain.
    fftBackward(hyperWork);
    for(i = 0; i < spatialCount; i++)
    {
        spat[i] *= hyper->spatial[i];
    }
    fftForward(hyperWork);
    for(i = 0; i < spectralCount; i++)
    {
        if(add)
            out[i] += spect[i];
        else
            out[i] = spect[i];
    }

    trace("Finished Hyper Diffusion\n");
}

/*
 * This method does not work. It is means to zero out the boundaries when the
 * a domain shift moves data through the periodic boundaries, so that fluid can 
 * rise through an infinite domain, but either the boundaries do not get
 * sanitized or else numerical instabilities ruin everything.
 */
void killBoundaries(PRECISION * in, complex<PRECISION> * out, int add, PRECISION factor)
{
    trace("Starting Boundary Forcing\n");
    
    int i,j,x,y,z;
    int index = 0;
    
    complex<PRECISION> * spect = hyperWork->spectral;
    PRECISION * spat = hyperWork->spatial;

    for(i = 0; i < spatialCount; i++)
    {
        spat[i] = -hyper->spatial[i]*factor*in[i];;
    }
    fftForward(hyperWork);
    for(i = 0; i < spectralCount; i++)
    {
        if(add)
            out[i] += spect[i];
        else
            out[i] = spect[i];
    }

    trace("Finished Boundary Forcing\n");
}

/*
 * Standard curl operation where:
 * 
 * Curl(A) = < dy Az - dz Ay, dz Ax - dx Az, dx Ay - dy Ax >
 * 
 * All wave modes are independent, so in and out can point to the same memory
 * as long as we are allowed to overwrite.
 */
void curl(p_vector in, p_vector out)
{
    int i,j,k;
    complex<PRECISION> dkx, dky, dkz;
    int index = 0;

    complex<PRECISION> * xin = in->x->spectral;
    complex<PRECISION> * yin = in->y->spectral;
    complex<PRECISION> * zin = in->z->spectral;

    complex<PRECISION> * xout = out->x->spectral;
    complex<PRECISION> * yout = out->y->spectral;
    complex<PRECISION> * zout = out->z->spectral;

    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);
                xout[index] = dky * zin[index] - dkz * yin[index];
                yout[index] = dkz * xin[index] - dkx * zin[index];
                zout[index] = dkx * yin[index] - dky * xin[index];

                index++;
            }
        }
    }
}

/*
 * Standard gradient of a scalar field where:
 * 
 * grad(A) = < dx A, dy A, dz A >
 * 
 * Again in and out can point to the same memory locations.
 */
void gradient(p_field in, p_vector out)
{
    int i,j,k;
    complex<PRECISION> dkx, dky, dkz;

    complex<PRECISION> * pin = in->spectral;
    complex<PRECISION> * outx = out->x->spectral;
    complex<PRECISION> * outy = out->y->spectral;
    complex<PRECISION> * outz = out->z->spectral;

    int index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);

                outx[index] = dkx * pin[index];
                outy[index] = dky * pin[index];
                outz[index] = dkz * pin[index];

                index++;
            }
        }
    }
}

/*
 * A dot B = Ax Bx + Ay By + Az Bz
 */
void dotProduct(p_vector one, p_vector two, p_field out)
{
    int i;

    PRECISION * onex = one->x->spatial;
    PRECISION * oney = one->y->spatial;
    PRECISION * onez = one->z->spatial;
    PRECISION * twox = two->x->spatial;
    PRECISION * twoy = two->y->spatial;
    PRECISION * twoz = two->z->spatial;
    PRECISION * pout = out->spatial;

    for(i = 0; i < spatialCount; i++)
    {
        pout[i] = onex[i] * twox[i];
        pout[i] += oney[i] * twoy[i];
        pout[i] += onez[i] * twoz[i];
    }
}

/*
 * A cross B = < Ay Bz - Az By, Az Bx - Ax Bz, Ax By - Ay Bx >
 */
void crossProduct(p_vector one, p_vector two, p_vector out)
{
    int i;

    PRECISION * onex = one->x->spatial;
    PRECISION * oney = one->y->spatial;
    PRECISION * onez = one->z->spatial;
    PRECISION * twox = two->x->spatial;
    PRECISION * twoy = two->y->spatial;
    PRECISION * twoz = two->z->spatial;
    PRECISION * outx = out->x->spatial;
    PRECISION * outy = out->y->spatial;
    PRECISION * outz = out->z->spatial;

    for(i = 0; i < spatialCount; i++)
    {
        outx[i] = oney[i] * twoz[i] - onez[i] * twoy[i];
        outy[i] = onez[i] * twox[i] - onex[i] * twoz[i];
        outz[i] = onex[i] * twoy[i] - oney[i] * twox[i];
    }
}

/*
 * div(a) = dx Ax + dy Ay + dz Az
 */
void divergence(p_vector in, p_field out)
{
    int i,j,k;
    complex<PRECISION> dkx,dky,dkz;
    int index = 0;

    complex<PRECISION> * o = out->spectral;
    complex<PRECISION> * x = in->x->spectral;
    complex<PRECISION> * y = in->y->spectral;
    complex<PRECISION> * z = in->z->spectral;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkz; k++)
            {
                dkz = dzFactor(k);

                o[index] = dkx * x[index];
                o[index] += dky * y[index];
                o[index] += dkz * z[index];

                index++;
            }
        }
    }

}

/*
 * Applied differentiation operation to a complex field, which reduces to
 * simply multiplying each element by its corresponding wavenyumber.
 * 
 * arithmetic = 0  : overwrite
 * arithmetic = 1  : +=
 * arithmetic = 2  : -=
 */
void partialX(complex<PRECISION> * in, complex<PRECISION> * out, int arithmetic)
{
    int i,j,k;
    int index = 0;
    complex<PRECISION> dk;

    if(arithmetic == 0)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            dk = dxFactor(i);
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    out[index] = dk*in[index];
                    index++;
                }
            }
        }
    }
    else if(arithmetic == 1)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            dk = dxFactor(i);
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    out[index] += dk*in[index];
                    index++;
                }
            }
        }
    }
    else if(arithmetic == 2)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            dk = dxFactor(i);
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    out[index] -= dk*in[index];
                    index++;
                }
            }
        }
    }
    else
    {
        error("Invalid arithmetic argument %d!!\n", arithmetic);
    }
}

/*
 * Applied differentiation operation to a complex field, which reduces to
 * simply multiplying each element by its corresponding wavenyumber.
 * 
 * arithmetic = 0  : overwrite
 * arithmetic = 1  : +=
 * arithmetic = 2  : -=
 */
void partialY(complex<PRECISION> * in, complex<PRECISION> * out, int arithmetic)
{
    int i,j,k;
    int index = 0;
    complex<PRECISION> dk;

    if(arithmetic == 0)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                dk = dyFactor(j);
                for(k = 0; k < ndkz; k++)
                {
                    out[index] = dk*in[index];
                    index++;
                }
            }
        }
    }
    else if(arithmetic == 1)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                dk = dyFactor(j);
                for(k = 0; k < ndkz; k++)
                {
                    out[index] += dk*in[index];
                    index++;
                }
            }
        }
    }
    else if(arithmetic == 2)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                dk = dyFactor(j);
                for(k = 0; k < ndkz; k++)
                {
                    out[index] -= dk*in[index];
                    index++;
                }
            }
        }
    }
    else
    {
        error("Invalid arithmetic argument %d!!\n", arithmetic);
    }
}

/*
 * Applied differentiation operation to a complex field, which reduces to
 * simply multiplying each element by its corresponding wavenyumber.
 * 
 * arithmetic = 0  : overwrite
 * arithmetic = 1  : +=
 * arithmetic = 2  : -=
 */
void partialZ(complex<PRECISION> * in, complex<PRECISION> * out, int arithmetic)
{
    int i,j,k;
    int index = 0;
    complex<PRECISION> dk;

    if(arithmetic == 0)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    dk = dzFactor(k);
                    out[index] = dk*in[index];
                    index++;
                }
            }
        }
    }
    else if(arithmetic == 1)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    dk = dzFactor(k);
                    out[index] += dk*in[index];
                    index++;
                }
            }
        }
    }
    else if(arithmetic == 2)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    dk = dzFactor(k);
                    out[index] -= dk*in[index];
                    index++;
                }
            }
        }
    }
    else
    {
        error("Invalid arithmetic argument %d!!\n", arithmetic);
    }
}

/*
 * Basic multiplication.  Can only be applied to fields in spatial coordinates
 * (no wave modes!)
 */
void multiply(PRECISION * one, PRECISION * two, PRECISION * out)
{
    int i;
    for(i = 0; i < spatialCount; i++)
    {
        out[i] = one[i] * two[i];
    }
}

/*
 * Basic addition routine for two complex fields.
 */
void plusEq(complex<PRECISION> * one, complex<PRECISION> * two)
{
    int i;
    for(i = 0; i < spectralCount; i++)
    {
        one[i] += two[i];
    }
}

/*
 * Basic subtraction routine for two complex fields.
 */
void minusEq(complex<PRECISION> * one, complex<PRECISION> * two)
{

    int i;
    for(i = 0; i < spectralCount; i++)
    {
        one[i] -= two[i];
    }
}

#include "Logs/LogInfo.h"
/*
 * To understand the following three routines, one must understand the data 
 * layout from an FFT operation.  For a basic 1D FFT, the first index of the 
 * resultant will be the mean, and each element afterwards can be thought of as
 * one wavenumber higher.  
 * 
 * However, due to aliasing there is an ambiguity, so if
 * there are N elements into the vector, index 7 can be thought of as wavenumber
 * 7, N+7, 2N+7, etc.  Since it makes sense from a physical standpoint to use
 * the smallest wavenumber modes, this means we will use index 7 as wavenumber
 * 7, but also means we will use index N-7 as wavenumber -7.  
 * 
 * The consequence of this is that our highest wavenumbers are actually 
 * stored in the center of the array.  So we always have to figure out which
 * half of the array we are in, and our current wavenumber is really the 
 * distance to the closest edge, not the distance to index 0.
 */
complex<PRECISION> dxFactor(int i)
{
    int k = i + my_kx->min;
    if(k >= dealias_kx.min)
    {
       k += 1 - 2 * dealias_kx.min;
    }

    complex<PRECISION> ret(0, 2 * PI * k / xmx);
    trace("kx index %d was found to be wave mode %d, %g\n", i, k, ret.imag() );

    return ret;
}

complex<PRECISION> dyFactor(int i)
{
    int k = i + my_ky->min;
    complex<PRECISION> ret(0, 2 * PI * k / ymx);
    trace("ky index %d was found to be wave mode %d, %g\n", i, k, ret.imag() );

    return ret;
}

complex<PRECISION> dzFactor(int i)
{
    int k = i;
    if(k >= dealias_kz.min)
    {
       k += 1 - 2 * dealias_kz.min;
    }
    complex<PRECISION> ret(0, 2 * PI * k / zmx);
    trace("kz index %d was found to be wave mode %d, %g\n", i, k, ret.imag());

    return ret;
}

/*
 * This routine shifts a given field f by a given displacement d.  This is done
 * by changing the phase of each wavemode appropriately.  This method works, but
 * is somewhat useless without the boundary sanitization routines that do not
 * currently work.
 */
void shiftField(displacement d, complex<PRECISION> * f)
{
    int i,j,k;
    complex<PRECISION> *edkx, *edky, *edkz;
    
    edkx = (complex<PRECISION>*)malloc(my_kx->width*sizeof(complex<PRECISION>));
    for(i=0; i < my_kx->width; i++)
        edkx[i] = exp(dxFactor(i)*d.dx);
    edky = (complex<PRECISION>*)malloc(my_ky->width*sizeof(complex<PRECISION>));
    for(i=0; i < my_ky->width; i++)
        edky[i] = exp(dyFactor(i)*d.dy);
    edkz = (complex<PRECISION>*)malloc(ndkz*sizeof(complex<PRECISION>));
    for(i=0; i < ndkz; i++)
        edkz[i] = exp(dzFactor(i)*d.dz);
    
    int index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < ndkz; k++)
            {
                f[index] *= edkx[i]*edky[j]*edkz[k];
                index++;
            }
        }
    }
    
    free(edkx);
    free(edky);
    free(edkz);
}

void shiftAvg(displacement d, complex<PRECISION> * f)
{
    int k;
    complex<PRECISION> dkz;
    
    int index = 0;
    for(k = 0; k < ndkz; k++)
    {
        dkz = dzFactor(k);

        f[index] *= exp(dkz*d.dz);

        index++;
    }
}


