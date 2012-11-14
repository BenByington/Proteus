#include "Physics.h"
#include "Numerics.h"
#include "State.h"
#include "Log.h"
#include "Communication.h"
#include "Environment.h"
#include "Field.h"
#include "TimeFunctions.h"

#include <string.h>
#include <math.h>

p_vector rhs = 0;
p_vector temp1 = 0;
p_vector temp2 = 0;

void calcForces();
void calcMomentum();
void calcTemp();
void calcMag();
void calcNewTimestep();
void step();

void eulerStep();
void AB2Step();
void AB3Step();

void iterate()
{
    calcNewTimestep();
    calcForces();
    step();

    //make sure our u,v,w terms are up to date, both spectral and spatial
    if(momEquation)
    {
        recomposeSolenoidal(u->sol, u->vec);
        fftBackward(u->vec->x);
        fftBackward(u->vec->y);
        fftBackward(u->vec->z);
    }

    if(magEquation)
    {
        recomposeSolenoidal(B->sol, B->vec);
        fftBackward(B->vec->x);
        fftBackward(B->vec->y);
        fftBackward(B->vec->z);
    }

    if(tEquation)
    {
        fftBackward(T);
    }
}

void calcNewTimestep()
{
    dt2 = dt1;
    dt1 = dt;
    dt = 0.01;  //just a trial.  This should be overwritten by at least one
                //of the following cases.

    if(momEquation && viscosity)
    {
        dt = safetyFactor * pow(fmin(fmin(dx,dy),dz),2) / Pr;
    }
    if(tEquation && tDiff)
    {
        PRECISION temp = safetyFactor * pow(fmin(fmin(dx,dy),dz),2);
        if(temp < dt)
            dt = temp;
    }
    if(magEquation && magDiff)
    {
        PRECISION temp = safetyFactor * pow(fmin(fmin(dx,dy),dz),2) * Pm / Pr;
        if(temp < dt)
            dt = temp;
    }
    if((momEquation && momAdvection) || (tEquation && tempAdvection) || (magEquation && magAdvect))
    {
        int i;
        //get maxV for stability condition.
        maxVel[0] = 0;
        maxVel[1] = 0;
        maxVel[2] = 0;
        PRECISION * x = u->vec->x->spatial;
        PRECISION * y = u->vec->y->spatial;
        PRECISION * z = u->vec->z->spatial;
        for(i = 0; i < spatialCount; i++)
        {
            if(fabs(x[i]) > maxVel[0])
                maxVel[0] = fabs(x[i]);
            if(fabs(y[i]) > maxVel[1])
                maxVel[1] = fabs(y[i]);
            if(fabs(z[i]) > maxVel[2])
                maxVel[2] = fabs(z[i]);
        }

        MPI_Allreduce(MPI_IN_PLACE, maxVel, 1, MPI_PRECISION, MPI_MAX, ccomm);
        MPI_Allreduce(MPI_IN_PLACE, maxVel+1, 1, MPI_PRECISION, MPI_MAX, ccomm);
        MPI_Allreduce(MPI_IN_PLACE, maxVel+2, 1, MPI_PRECISION, MPI_MAX, ccomm);

        trace("Max VeL %g %g %g\n", maxVel[0], maxVel[1], maxVel[2]);
        
        if(maxVel[0] != 0)
        {
            if(safetyFactor * dx / maxVel[0] < dt)
                dt = safetyFactor * dx / maxVel[0];
        }
        if(maxVel[1] != 0)
        {
            if(safetyFactor * dy / maxVel[1] < dt)
                dt = safetyFactor * dy / maxVel[1];
        }
        if(maxVel[2] != 0)
        {
            if(safetyFactor * dz / maxVel[2] < dt)
                dt = safetyFactor * dz / maxVel[2];
        }
    }
    trace("time step for iteration %d is %g\n",iteration, dt);

}

void calcForces()
{
    debug("Calculating forces\n",0);

    //cycle the force pointers
    complex PRECISION * temp;
    if(momEquation)
    {
        temp = u->sol->poloidal->force3;
        u->sol->poloidal->force3 = u->sol->poloidal->force2;
        u->sol->poloidal->force2 = u->sol->poloidal->force1;
        u->sol->poloidal->force1 = temp;

        temp = u->sol->toroidal->force3;
        u->sol->toroidal->force3 = u->sol->toroidal->force2;
        u->sol->toroidal->force2 = u->sol->toroidal->force1;
        u->sol->toroidal->force1 = temp;

        temp = u->sol->mean_xf3;
        u->sol->mean_xf3 = u->sol->mean_xf2;
        u->sol->mean_xf2 = u->sol->mean_xf1;
        u->sol->mean_xf1 = temp;

        temp = u->sol->mean_yf3;
        u->sol->mean_yf3 = u->sol->mean_yf2;
        u->sol->mean_yf2 = u->sol->mean_yf1;
        u->sol->mean_yf1 = temp;
    }

    if(magEquation)
    {
        temp = B->sol->poloidal->force3;
        B->sol->poloidal->force3 = B->sol->poloidal->force2;
        B->sol->poloidal->force2 = B->sol->poloidal->force1;
        B->sol->poloidal->force1 = temp;

        temp = B->sol->toroidal->force3;
        B->sol->toroidal->force3 = B->sol->toroidal->force2;
        B->sol->toroidal->force2 = B->sol->toroidal->force1;
        B->sol->toroidal->force1 = temp;

        temp = B->sol->mean_xf3;
        B->sol->mean_xf3 = B->sol->mean_xf2;
        B->sol->mean_xf2 = B->sol->mean_xf1;
        B->sol->mean_xf1 = temp;

        temp = B->sol->mean_yf3;
        B->sol->mean_yf3 = B->sol->mean_yf2;
        B->sol->mean_yf2 = B->sol->mean_yf1;
        B->sol->mean_yf1 = temp;
    }

    if(tEquation)
    {
        temp = T->force3;
        T->force3 = T->force2;
        T->force2 = T->force1;
        T->force1 = temp;
    }

    if(momEquation)
        calcMomentum();

    if(tEquation)
        calcTemp();

    if(magEquation)
        calcMag();

    debug("Forces done\n",0);
}

void calcMomentum()
{
    debug("Calculating momentum forces\n",0);

    int i,j,k;
    int index;

    if(viscosity)
    {
        laplacian(u->vec->x->spectral, rhs->x->spectral, 0, Pr);
        laplacian(u->vec->y->spectral, rhs->y->spectral, 0, Pr);
        laplacian(u->vec->z->spectral, rhs->z->spectral, 0, Pr);
    }
    else
    {
        //make sure we don't start with garbage
        memset(rhs->x->spectral, 0, spectralCount * sizeof(complex PRECISION));
        memset(rhs->y->spectral, 0, spectralCount * sizeof(complex PRECISION));
        memset(rhs->z->spectral, 0, spectralCount * sizeof(complex PRECISION));
    }

    //static forcing is currently only in the u direction and a function of y and z
    if(momStaticForcing)
    {
        complex PRECISION * xfield = rhs->x->spectral;
        complex PRECISION * ffield = forceField->spectral;
        index = 0;
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    xfield[index] += ffield[index];
                    index++;
                }
            }
        }
    }

    if(momTimeForcing)
    {
        fillTimeField(temp1, MOMENTUM);

        complex PRECISION * xfield = rhs->x->spectral;
        complex PRECISION * yfield = rhs->y->spectral;
        complex PRECISION * zfield = rhs->z->spectral;
        complex PRECISION * xforce = temp1->x->spectral;
        complex PRECISION * yforce = temp1->y->spectral;
        complex PRECISION * zforce = temp1->z->spectral;

        int i;
        for(i = 0; i < spectralCount; i++)
        {
            xfield[i] += xforce[i];
            yfield[i] += yforce[i];
            zfield[i] += zforce[i];
        }
    }

    
    if(momAdvection)
    {
        p_field tense = temp1->x;

        multiply(u->vec->x->spatial, u->vec->x->spatial, tense->spatial);
        fftForward(tense);
        partialX(tense->spectral, rhs->x->spectral, 2);

        multiply(u->vec->y->spatial, u->vec->y->spatial, tense->spatial);
        fftForward(tense);
        partialY(tense->spectral, rhs->y->spectral, 2);

        multiply(u->vec->z->spatial, u->vec->z->spatial, tense->spatial);
        fftForward(tense);
        partialZ(tense->spectral, rhs->z->spectral, 2);

        multiply(u->vec->x->spatial, u->vec->y->spatial, tense->spatial);
        fftForward(tense);
        partialY(tense->spectral, rhs->x->spectral, 2);
        partialX(tense->spectral, rhs->y->spectral, 2);

        multiply(u->vec->x->spatial, u->vec->z->spatial, tense->spatial);
        fftForward(tense);
        partialZ(tense->spectral, rhs->x->spectral, 2);
        partialX(tense->spectral, rhs->z->spectral, 2);

        multiply(u->vec->y->spatial, u->vec->z->spatial, tense->spatial);
        fftForward(tense);
        partialZ(tense->spectral, rhs->y->spectral, 2);
        partialY(tense->spectral, rhs->z->spectral, 2);      
    }


    if(lorentz)
    {
        p_field tense = temp1->x;
        p_vector lor = temp2;

        multiply(B->vec->x->spatial, B->vec->x->spatial, tense->spatial);
        fftForward(tense);
        partialX(tense->spectral, lor->x->spectral, 0);

        multiply(B->vec->y->spatial, B->vec->y->spatial, tense->spatial);
        fftForward(tense);
        partialY(tense->spectral, lor->y->spectral, 0);

        multiply(B->vec->z->spatial, B->vec->z->spatial, tense->spatial);
        fftForward(tense);
        partialZ(tense->spectral, lor->z->spectral, 0);

        multiply(B->vec->x->spatial, B->vec->y->spatial, tense->spatial);
        fftForward(tense);
        partialY(tense->spectral, lor->x->spectral, 1);
        partialX(tense->spectral, lor->y->spectral, 1);

        multiply(B->vec->x->spatial, B->vec->z->spatial, tense->spatial);
        fftForward(tense);
        partialZ(tense->spectral, lor->x->spectral, 1);
        partialX(tense->spectral, lor->z->spectral, 1);

        multiply(B->vec->y->spatial, B->vec->z->spatial, tense->spatial);
        fftForward(tense);
        partialZ(tense->spectral, lor->y->spectral, 1);
        partialY(tense->spectral, lor->z->spectral, 1);

        //TODO Find a routine to change to include this factor
        complex PRECISION * px = lor->x->spectral;
        complex PRECISION * py = lor->y->spectral;
        complex PRECISION * pz = lor->z->spectral;
        int i;
        for(i = 0; i < spectralCount; i++)
        {
            px[i] *= alpha;
            py[i] *= alpha;
            pz[i] *= alpha;
        }

        plusEq(rhs->x->spectral, px);
        plusEq(rhs->y->spectral, py);
        plusEq(rhs->z->spectral, pz);
    }

    if(buoyancy)
    {
        complex PRECISION * zfield = rhs->z->spectral;
        complex PRECISION * tfield = T->spectral;

        PRECISION factor =  Ra * Pr;
        index = 0;
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    zfield[index] += factor * tfield[index];

                    index++;
                }
            }
        }
    }
    
    //curl it so we can avoid dealing with the pressure term
    curl(rhs, temp1);

    decomposeCurlSolenoidal(u->sol, temp1, 1);
    debug("Momentum forces done\n",0);
}

void calcMag()
{
    int i,j,k;
    int index;
    debug("Calculating Magnetic forces\n",0);

    if(magDiff)
    {
        laplacian(B->vec->x->spectral, rhs->x->spectral, 0, Pr/Pm);
        laplacian(B->vec->y->spectral, rhs->y->spectral, 0, Pr/Pm);
        laplacian(B->vec->z->spectral, rhs->z->spectral, 0, Pr/Pm);
    }
    else
    {
        //make sure we don't start with garbage
        memset(rhs->x->spectral, 0, spectralCount * sizeof(complex PRECISION));
        memset(rhs->y->spectral, 0, spectralCount * sizeof(complex PRECISION));
        memset(rhs->z->spectral, 0, spectralCount * sizeof(complex PRECISION));
    }

    //static forcing is currently only in the x direction and a function of y and z
    if(magStaticForcing)
    {
        complex PRECISION * xfield = rhs->x->spectral;
        complex PRECISION * ffield = magForceField->spectral;
        index = 0;
        for(i = 0; i < my_kx->width; i++)
        {
            for(j = 0; j < my_ky->width; j++)
            {
                for(k = 0; k < ndkz; k++)
                {
                    xfield[index] += ffield[index];
                    index++;
                }
            }
        }
    }
    
    if(kinematic)
    {
        fillTimeField(u->vec, KINEMATIC);
    }

    if(magAdvect)
    {
        p_vector uxb = temp1;
        p_vector cuxb = temp2;

        crossProduct(u->vec, B->vec, uxb);
        fftForward(uxb->x);
        fftForward(uxb->y);
        fftForward(uxb->z);

        curl(uxb, cuxb);

        plusEq(rhs->x->spectral, cuxb->x->spectral);
        plusEq(rhs->y->spectral, cuxb->y->spectral);
        plusEq(rhs->z->spectral, cuxb->z->spectral);
    }

    if(magTimeForcing)
    {
        fillTimeField(temp1, MAGNETIC);

        plusEq(rhs->x->spectral, temp1->x->spectral);
        plusEq(rhs->y->spectral, temp1->y->spectral);
        plusEq(rhs->z->spectral, temp1->z->spectral);
    }

    decomposeSolenoidal(B->sol, rhs, 1);
    debug("Magnetic forces done\n",0);
}

void calcTemp()
{
    complex PRECISION * forces = T->force1;
    if(tDiff)
    {
        laplacian(T->spectral, forces, 0, 1.0);
    }
    else
    {
        memset(forces, 0, spectralCount * sizeof(complex PRECISION));
    }

    if(tempAdvection)
    {
        //advect the background profile
        plusEq(forces, u->vec->z->spectral);

        //advect the perturbations
        p_vector flux = temp1;
        multiply(u->vec->x->spatial, T->spatial, flux->x->spatial);
        multiply(u->vec->y->spatial, T->spatial, flux->y->spatial);
        multiply(u->vec->z->spatial, T->spatial, flux->z->spatial);

        fftForward(flux->x);
        fftForward(flux->y);
        fftForward(flux->z);

        p_field advect = temp2->x;
        divergence(flux, advect);

        minusEq(forces, advect->spectral);

    }
}

void initPhysics()
{
    temp1 = newVector(SPEC | SPAT);
    temp2 = newVector(SPEC | SPAT);
    rhs = newVector(SPEC);
}

void finalizePhysics()
{
    deleteVector(&temp1);
    deleteVector(&temp2);
    deleteVector(&rhs);
}

void step()
{
    if(iteration == 1)
    {
        eulerStep();
    }
    else if(iteration == 2)
    {
        AB2Step();
    }
    else
    {
        AB3Step();
    }
    elapsedTime += dt;
}

void eulerStep()
{
    int i;
    
    complex PRECISION * func;
    complex PRECISION * f1;

    if(momEquation)
    {
        func = u->sol->poloidal->spectral;
        f1 = u->sol->poloidal->force1;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += dt * f1[i];
        }

        func = u->sol->toroidal->spectral;
        f1 = u->sol->toroidal->force1;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += dt * f1[i];
        }

        func = u->sol->mean_x;
        f1 = u->sol->mean_xf1;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += dt * f1[i];
        }

        func = u->sol->mean_y;
        f1 = u->sol->mean_yf1;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += dt * f1[i];
        }
    }

    if(magEquation)
    {
        func = B->sol->poloidal->spectral;
        f1 = B->sol->poloidal->force1;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += dt * f1[i];
        }

        func = B->sol->toroidal->spectral;
        f1 = B->sol->toroidal->force1;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += dt * f1[i];
        }

        func = B->sol->mean_x;
        f1 = B->sol->mean_xf1;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += dt * f1[i];
        }

        func = B->sol->mean_y;
        f1 = B->sol->mean_yf1;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += dt * f1[i];
        }
    }

    if(tEquation)
    {
        func = T->spectral;
        f1 = T->force1;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += dt * f1[i];
        }
    }
}

void AB2Step()
{
    int i;

    complex PRECISION * func;
    complex PRECISION * f1;
    complex PRECISION * f2;

    PRECISION c0 = dt * (0.5 * dt / dt1 + 1);
    PRECISION c1 = -0.5 * dt * dt / dt1;

    if(momEquation)
    {
        func = u->sol->poloidal->spectral;
        f1 = u->sol->poloidal->force1;
        f2 = u->sol->poloidal->force2;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }

        func = u->sol->toroidal->spectral;
        f1 = u->sol->toroidal->force1;
        f2 = u->sol->toroidal->force2;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }

        func = u->sol->mean_x;
        f1 = u->sol->mean_xf1;
        f2 = u->sol->mean_xf2;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }

        func = u->sol->mean_y;
        f1 = u->sol->mean_yf1;
        f2 = u->sol->mean_yf2;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }
    }

    if(magEquation)
    {
        func = B->sol->poloidal->spectral;
        f1 = B->sol->poloidal->force1;
        f2 = B->sol->poloidal->force2;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }

        func = B->sol->toroidal->spectral;
        f1 = B->sol->toroidal->force1;
        f2 = B->sol->toroidal->force2;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }

        func = B->sol->mean_x;
        f1 = B->sol->mean_xf1;
        f2 = B->sol->mean_xf2;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }

        func = B->sol->mean_y;
        f1 = B->sol->mean_yf1;
        f2 = B->sol->mean_yf2;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }
    }

    if(tEquation)
    {
        func = T->spectral;
        f1 = T->force1;
        f2 = T->force2;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i];
        }
    }
}

void AB3Step()
{
    int i;

    complex PRECISION * func;
    complex PRECISION * f1;
    complex PRECISION * f2;
    complex PRECISION * f3;

    PRECISION c0 =  dt + (dt/dt1)*(dt/(dt1+dt2))*(dt/3.0 + 0.5*(2*dt1+ dt2));
    PRECISION c1 = -(dt/dt1)*(dt/(dt2))*(dt/3.0 + 0.5*(dt1+dt2));
    PRECISION c2 = (dt/(dt1 + dt2))*(dt/(dt2))*(dt/3.0 + 0.5*dt1);

    if(momEquation)
    {
        func = u->sol->poloidal->spectral;
        f1 = u->sol->poloidal->force1;
        f2 = u->sol->poloidal->force2;
        f3 = u->sol->poloidal->force3;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }

        func = u->sol->toroidal->spectral;
        f1 = u->sol->toroidal->force1;
        f2 = u->sol->toroidal->force2;
        f3 = u->sol->toroidal->force3;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }

        func = u->sol->mean_x;
        f1 = u->sol->mean_xf1;
        f2 = u->sol->mean_xf2;
        f3 = u->sol->mean_xf3;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }

        func = u->sol->mean_y;
        f1 = u->sol->mean_yf1;
        f2 = u->sol->mean_yf2;
        f3 = u->sol->mean_yf3;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }
    }

    if(magEquation)
    {
        func = B->sol->poloidal->spectral;
        f1 = B->sol->poloidal->force1;
        f2 = B->sol->poloidal->force2;
        f3 = B->sol->poloidal->force3;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }

        func = B->sol->toroidal->spectral;
        f1 = B->sol->toroidal->force1;
        f2 = B->sol->toroidal->force2;
        f3 = B->sol->toroidal->force3;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }

        func = B->sol->mean_x;
        f1 = B->sol->mean_xf1;
        f2 = B->sol->mean_xf2;
        f3 = B->sol->mean_xf3;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }

        func = B->sol->mean_y;
        f1 = B->sol->mean_yf1;
        f2 = B->sol->mean_yf2;
        f3 = B->sol->mean_yf3;
        for(i = 0; i < ndkz; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }
    }

    if(tEquation)
    {
        func = T->spectral;
        f1 = T->force1;
        f2 = T->force2;
        f3 = T->force3;
        for(i = 0; i < spectralCount; i++)
        {
            func[i] += c0 * f1[i] + c1 * f2[i] + c2 * f3[i];
        }
    }
}

