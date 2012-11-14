#include "Numerics.h"
#include "Log.h"
#include "Environment.h"
#include "Communication.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

void testPT()
{
    int i,j,k,l;
    int index = 0;
    double ampM = 3;
    complex double dkx,dky,dkz;
    p_solenoid s = newSolenoid();
    p_vector v = newVector(SPEC);
    p_vector v2 = newVector(SPEC);

    //set the z field
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
                    v->z->spectral[index] = ampM * (((double)rand()) / RAND_MAX + I * ((double)rand()) / RAND_MAX);
                index++;
            }
        }
    }
    
    //set the y field
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
                    v->y->spectral[index] = ampM * (((double)rand()) / RAND_MAX + I * ((double)rand()) / RAND_MAX);
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
                    v->x->spectral[index] = ampM * (((double)rand()) / RAND_MAX + I * ((double)rand()) / RAND_MAX);
                }
                //everything else must cancle with other fields
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


    decomposeSolenoidal(s, v, 0);
    recomposeSolenoidal(s, v2);

    index = 0;
    for(i = 0; i < my_kx->width; i++)
    {
        for(j = 0; j < my_ky->width; j++)
        {
            for(k = 0; k < ndkz; k++)
            {
                double diff = cabs(v2->z->spectral[index] - v->z->spectral[index]);
                if(diff > 1e-10)
                    fprintf(stderr,"diffz: %d %d %d %d %g\n", grank, i + my_kx->min, j + my_ky->min, k, diff);

                diff = cabs(v2->y->spectral[index] - v->y->spectral[index]);
                if(diff > 1e-10)
                    fprintf(stderr,"diffy: %d %d %d %d %g\n", grank, i + my_kx->min, j + my_ky->min, k, diff);

                diff = cabs(v2->x->spectral[index] - v->x->spectral[index]);
                if(diff > 1e-10)
                    fprintf(stderr,"diffx: %d %d %d %d %g\n", grank, i + my_kx->min, j + my_ky->min, k, diff);
                index++;
            }
        }
    }

}

void decomposeSolenoidal(p_solenoid s, p_vector v, int force)
{
    //TODO: do some code duplication that will allow us to avoid checking
    //for 0 wave modes every iteration of the loops
    int i,j,k;
    debug("Beginning decomposition of solenoidal field\n",0);

    if(crank == 0)
    {
        s->mean_z = v->z->spectral[0];
    }

    debug("Calculating Poloidal field\n",0);

    complex double dkx,dky,dkz;
    complex double * pz = v->z->spectral;
    complex double * py = v->y->spectral;
    complex double * px = v->x->spectral;
    complex double * ppol;
    complex double * ptor;
    complex double * xmean;
    complex double * ymean;

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

    debug("Calculating Toroidal Field\n",0);
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

void decomposeCurlSolenoidal(p_solenoid s, p_vector v, int force)
{
    //TODO: do some code duplication that will allow us to avoid checking
    //for 0 wave modes every iteration of the loops
    int i,j,k;
    debug("Beginning decomposition of curled solenoidal field\n",0);

    if(crank == 0)
    {
        s->mean_z = 0;
    }

    debug("Calculating Toroidal field\n",0);

    complex double dkx,dky,dkz;
    complex double * pz = v->z->spectral;
    complex double * py = v->y->spectral;
    complex double * px = v->x->spectral;
    complex double * ppol;
    complex double * ptor;
    complex double * xmean;
    complex double * ymean;

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

    debug("Calculating Poloidal Field\n",0);
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

void recomposeSolenoidal(p_solenoid s, p_vector v)
{
    int i,j,k;
    complex double dkx,dky,dkz;
    complex double * pz = v->z->spectral;
    complex double * py = v->y->spectral;
    complex double * px = v->x->spectral;
    complex double * ppol = s->poloidal->spectral;
    complex double * ptor = s->toroidal->spectral;
    complex double * xmean = s->mean_x;
    complex double * ymean = s->mean_y;

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

//TODO: Verify that it is safe to have both in and out reference the same
//structure when in wont be needed any longer
extern void laplacian(complex double * in, complex double * out, int add, double factor)
{
    trace("Starting Laplacian\n",0);
    
    int i,j,k;
    int index = 0;
    complex double dkx,dky,dkz;

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

    trace("Finished Laplacian\n",0);
}

extern void curl(p_vector in, p_vector out)
{
    int i,j,k;
    complex double dkx, dky, dkz;
    int index = 0;

    complex double * xin = in->x->spectral;
    complex double * yin = in->y->spectral;
    complex double * zin = in->z->spectral;

    complex double * xout = out->x->spectral;
    complex double * yout = out->y->spectral;
    complex double * zout = out->z->spectral;

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

extern void gradient(p_field in, p_vector out)
{
    int i,j,k;
    complex double dkx, dky, dkz;

    complex double * pin = in->spectral;
    complex double * outx = out->x->spectral;
    complex double * outy = out->y->spectral;
    complex double * outz = out->z->spectral;

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

extern void dotProduct(p_vector one, p_vector two, p_field out)
{
    int i;

    double * onex = one->x->spatial;
    double * oney = one->y->spatial;
    double * onez = one->z->spatial;
    double * twox = two->x->spatial;
    double * twoy = two->y->spatial;
    double * twoz = two->z->spatial;
    double * pout = out->spatial;

    for(i = 0; i < spatialCount; i++)
    {
        pout[i] = onex[i] * twox[i];
        pout[i] += oney[i] * twoy[i];
        pout[i] += onez[i] * twoz[i];
    }
}

extern void crossProduct(p_vector one, p_vector two, p_vector out)
{
    int i;

    double * onex = one->x->spatial;
    double * oney = one->y->spatial;
    double * onez = one->z->spatial;
    double * twox = two->x->spatial;
    double * twoy = two->y->spatial;
    double * twoz = two->z->spatial;
    double * outx = out->x->spatial;
    double * outy = out->y->spatial;
    double * outz = out->z->spatial;

    for(i = 0; i < spatialCount; i++)
    {
        outx[i] = oney[i] * twoz[i] - onez[i] * twoy[i];
        outy[i] = onez[i] * twox[i] - onex[i] * twoz[i];
        outz[i] = onex[i] * twoy[i] - oney[i] * twox[i];
    }
}

extern void divergence(p_vector in, p_field out)
{
    int i,j,k;
    complex double dkx,dky,dkz;
    int index = 0;

    complex double * o = out->spectral;
    complex double * x = in->x->spectral;
    complex double * y = in->y->spectral;
    complex double * z = in->z->spectral;
    for(i = 0; i < my_kx->width; i++)
    {
        dkx = dxFactor(i);
        for(j = 0; j < my_ky->width; j++)
        {
            dky = dyFactor(j);
            for(k = 0; k < ndkx; k++)
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

extern void partialX(complex double * in, complex double * out, int arithmetic)
{
    int i,j,k;
    int index = 0;
    complex double dk;

    if(arithmetic == 0)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            dk = dxFactor(i);
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    out[i] = dk*in[i];
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
                    out[i] += dk*in[i];
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
                    out[i] -= dk*in[i];
                    index++;
                }
            }
        }
    }
    else
        error("Invalid arithmetic argument %d!!\n", arithmetic);
}

extern void partialY(complex double * in, complex double * out, int arithmetic)
{
    int i,j,k;
    int index = 0;
    complex double dk;

    if(arithmetic == 0)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                dk = dyFactor(j);
                for(k = 0; k < ndkz; k++)
                {
                    out[i] = dk*in[i];
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
                    out[i] += dk*in[i];
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
                    out[i] -= dk*in[i];
                    index++;
                }
            }
        }
    }
    else
        error("Invalid arithmetic argument %d!!\n", arithmetic);

}

extern void partialZ(complex double * in, complex double * out, int arithmetic)
{
    int i,j,k;
    int index = 0;
    complex double dk;

    if(arithmetic == 0)
    {
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    dk = dzFactor(k);
                    out[i] = dk*in[i];
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
                    out[i] += dk*in[i];
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
                    out[i] -= dk*in[i];
                    index++;
                }
            }
        }
    }
    else
        error("Invalid arithmetic argument %d!!\n", arithmetic);
}

extern void multiply(double * one, double * two, double * out)
{
    int i;
    for(i = 0; i < spatialCount; i++)
    {
        out[i] = one[i] * two[i];
    }
}

extern void plusEq(complex double * one, complex double * two)
{

    int i;
    for(i = 0; i < spectralCount; i++)
    {
        one[i] += two[i];
    }
}

extern void minusEq(complex double * one, complex double * two)
{

    int i;
    for(i = 0; i < spectralCount; i++)
    {
        one[i] -= two[i];
    }
}

#include "LogInfo.h"
extern complex double dxFactor(int i)
{
    int k = i + my_kx->min;
    if(k >= dealias_kx.min)
    {
       k += 1 - 2 * dealias_kx.min;
    }

    trace("kx index %d was found to be wave mode %d, %g\n", i, k, __imag__ (I * 2 * PI * k / xmx) );

    return I * 2 * PI * k / xmx;
}

extern complex double dyFactor(int i)
{
    int k = i + my_ky->min;
    complex double ret = I * 2 * PI * k / ymx;
    trace("ky index %d was found to be wave mode %d, %g\n", i, k, __imag__ ret );

    return ret;
}

extern complex double dzFactor(int i)
{
    int k = i;
    if(k >= dealias_kz.min)
    {
       k += 1 - 2 * dealias_kz.min;
    }
    trace("kz index %d was found to be wave mode %d, %g\n", i, k, __imag__ (I * 2 * PI * k / xmx) );

    return I * 2 * PI * k / zmx;
}

