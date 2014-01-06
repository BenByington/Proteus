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

#include "Physics.h"
#include "Numerics.h"
#include "State.h"
#include "Logs/Log.h"
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

/* 
 * This is one of the few methods available externally.  Here we simply 
 * calculate the timestep to use for this iteration, calculate the new batch
 * of forces, and then propagate our variables forward through time.
 */
void iterate()
{
    calcNewTimestep();
    calcForces();
    step();

    /*
     * This is an experimental and only partially functional attempt to recenter
     * the domain upon a quantity of interest.  
     */
    if(recentering != NOCENTERING)
    {
        displacement d;
        
        if(recentering == BYMAXCENTER)
            d = displacementByCenter();
        
        if(momEquation)
        {
            shiftField(d, u->sol->toroidal->spectral);
            shiftField(d, u->sol->toroidal->force1);
            shiftField(d, u->sol->toroidal->force2);
            
            shiftField(d, u->sol->poloidal->spectral);
            shiftField(d, u->sol->poloidal->force1);
            shiftField(d, u->sol->poloidal->force2);
            
            shiftAvg(d, u->sol->mean_x);
            shiftAvg(d, u->sol->mean_xf1);
            shiftAvg(d, u->sol->mean_xf2);
            
            shiftAvg(d, u->sol->mean_y);
            shiftAvg(d, u->sol->mean_yf1);
            shiftAvg(d, u->sol->mean_yf2);
        }
        if(magEquation)
        {
            shiftField(d, B->sol->toroidal->spectral);
            shiftField(d, B->sol->toroidal->force1);
            shiftField(d, B->sol->toroidal->force2);
            
            shiftField(d, B->sol->poloidal->spectral);
            shiftField(d, B->sol->poloidal->force1);
            shiftField(d, B->sol->poloidal->force2);
            
            shiftAvg(d, B->sol->mean_x);
            shiftAvg(d, B->sol->mean_xf1);
            shiftAvg(d, B->sol->mean_xf2);
            
            shiftAvg(d, B->sol->mean_y);
            shiftAvg(d, B->sol->mean_yf1);
            shiftAvg(d, B->sol->mean_yf2);
        }
        if(tEquation)
        {
            shiftField(d, T->spectral);
            shiftField(d, T->force1);
            shiftField(d, T->force2);
        }
    }
    
    //make sure our state variables are up to date, both spectral and spatial
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

/*
 * This routine calculates the dt to use for the next time step.  This dt is
 * intentionally variable so that we can take as large of a time step as is
 * stable, depending on the state of the system.  There are two types of
 * guesses included in the routine below, both of which essentially boil down
 * to making sure that no information travels further than one grid cell during
 * an iteration.  To take larger steps one would need implicit methods.
 * 
 * 1.   The diffusion timescale in the x direction with a diffusion coefficient
 *      of D is dx**2 / D.  For each diffusion term in the problem, the most 
 *      restrictive dt is chosen.
 * 
 * 2.   For advection terms, we simply ensure that when using the peak 
 *      velocities in the problem, a passive tracer particle cannot move from
 *      one boundary to the next.
 * 
 * While there are sets of equations where these guidelines are actual stability
 * criterion, in this more complicated set of equations they are simply good
 * ideas.  As such, there is also a safety factor multiplied by all of our
 * estimates in order to help maintain stability.  In the past, a safety factor
 * of rougly 0.02 has been used, though it may be possible to safely raise this
 * further.
 */
void calcNewTimestep()
{
    //Shuffle along our old dt's.  We need to record what the past two were
    //for our AB3 routine.
    dt2 = dt1;
    dt1 = dt;
    dt = 0.01;  //just a trial.  This should be overwritten by at least one
                //of the following cases.

    PRECISION min_d2 = pow(fmin(fmin(dx,dy),dz),2);
    
    if(momEquation && viscosity)
    {
        dt = safetyFactor * min_d2 / Pr;
    }
    if(tEquation && tDiff)
    {
        PRECISION temp = safetyFactor * min_d2;
        if(temp < dt)
            dt = temp;
    }
    if(magEquation && magDiff)
    {
        PRECISION temp = safetyFactor * min_d2 * Pm / Pr;
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
        
        //Each direction may have a different restraint.  Take the most 
        //restrictive!
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

/*
 * This is just an entry point.  All this routine does is shuffle the force
 * fields, so that we can save the last two forcing evaluations, and call
 * the routines in charge of each individual equation.
 */
void calcForces()
{
    debug("Calculating forces\n");

    //cycle the force pointers for any active equation.
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

    //real force calculations are in these methods.
    if(momEquation)
        calcMomentum();

    if(tEquation)
        calcTemp();

    if(magEquation)
        calcMag();
   
    debug("Forces done\n");
}

/*
 * This routine is in charge of the momentum equation.  Virtually all
 * of the terms can be enabled or disabled by parameters read in through the
 * configuration file.  This equation has the form:
 * 
 * du/dt = div(alpha BB - uu - P) + (Ra T + B^2 / beta) hat z + Pr del^2 u + F_u 
 * 
 * Since the velocity field is incompressible, the pressure term doesn't matter
 * and is eliminated by taking the curl of the right-hand side.  This is then
 * decomposed into poloidal and toroidal components, which are then used in the
 * time integration.
 */
void calcMomentum()
{
    debug("Calculating momentum forces\n");

    int i,j,k;
    int index;

    if(viscosity)
    {
        //First argument is the field we take the laplacian of.
        //Second argument is where the result is stored.
        //The 0 in the third argument means we overwrite the destination array
        //Pr is the coefficient for this diffusion term.
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
    
    //Apply hyper diffusion to the boundaries
    //This does not currently work and is disabled by default!
    if(sanitize)
    {
        killBoundaries(u->vec->x->spatial, rhs->x->spectral, 1, 100*Pr);
        killBoundaries(u->vec->y->spatial, rhs->y->spectral, 1, 100*Pr);
        killBoundaries(u->vec->z->spatial, rhs->z->spectral, 1, 100*Pr);
    }

    //static forcing is read in from a file and currently only in the u 
    //direction as a function of y and z (to remove nonlinear advection)
    if(momStaticForcing)
    {
        complex PRECISION * xfield = rhs->x->spectral;
        complex PRECISION * ffield = forceField->spectral;
        index = 0;
        for(i = 0; i < spectralCount; i++)
        {
            xfield[i] += ffield[i];
        }
    }

    //
    if(momTimeForcing)
    {
        //Evaluate the force function at the current time.
        fillTimeField(temp1, MOMENTUM);

        plusEq(rhs->x->spectral, temp1->x->spectral);
        plusEq(rhs->y->spectral, temp1->y->spectral);
        plusEq(rhs->z->spectral, temp1->z->spectral);
    }

    
    if(momAdvection)
    {
        p_field tense = temp1->x;

	    //Note here, because I already forgot once.  a 2 as the third
	    //parameter makes things behave as a -= operation!
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

        //The third parameter as a 0 means we overwrite the destination array.
        //The third parameter as a 1 means it behaves as a += operation.
        
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

        //TODO: Find a routine to change to include this factor
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
    
    if(magBuoy)
    {
        p_field B2 = temp1->x;
        dotProduct(B->vec,B->vec,B2);
        fftForward(B2);
        for(i = 0; i < spectralCount; i++)
        {
            B2->spectral[i] *= magBuoyScale;
        }
        plusEq(rhs->z->spectral, B2->spectral);
    }

    if(buoyancy)
    {
        complex PRECISION * zfield = rhs->z->spectral;
        complex PRECISION * tfield = T->spectral;

        PRECISION factor =  Ra * Pr;
        index = 0;
        for(i = 0; i < spectralCount; i++)
        {
            zfield[i] += factor * tfield[i];
        }
    }
    
    //curl it so we can avoid dealing with the pressure term
    curl(rhs, temp1);

    //The third parameter as a 1 means we store the result in the force
    //arrays for u->sol.
    decomposeCurlSolenoidal(u->sol, temp1, 1);
    debug("Momentum forces done\n");
}

/*
 * This routine is in charge of the magnetic induction equation.  Virtually all
 * of the terms can be enabled or disabled by parameters read in through the
 * configuration file.  This equation has the form:
 * 
 * dB/dt = curl(u cross B) + Pm/Pr*del^2 B + F_B 
 * 
 * Again we have an incompressible field, so after evaluating the force in 
 * vector notation, we decompose it into poloidal and toroidal scalar fields for
 * the time integration.
 */
void calcMag()
{
    int i,j,k;
    int index;
    debug("Calculating Magnetic forces\n");

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
    
    //Apply hyper diffusion to the boundaries.  Again, this does not currently
    //work!
    if(sanitize)
    {
        killBoundaries(B->vec->x->spatial, rhs->x->spectral, 1, 100*Pr/Pm);
        killBoundaries(B->vec->y->spatial, rhs->y->spectral, 1, 100*Pr/Pm);
        killBoundaries(B->vec->z->spatial, rhs->z->spectral, 1, 100*Pr/Pm);
    }

    //static forcing is currently only in the x direction and a function of y and z
    if(magStaticForcing)
    {
        complex PRECISION * xfield = rhs->x->spectral;
        complex PRECISION * ffield = magForceField->spectral;
        index = 0;
        for(i = 0; i < spectralCount; i++)
        {
            xfield[i] += ffield[i];
        }
    }
    
    //If we are doing the kinematic problem then the velocity field is
    //specified ahead of time.  Note, this conflicts with the momentum equation
    //being enabled.  If both momentum equation and kinematic are enabled, then
    //we will still be doing all the work of both, but right here we will erase
    //any work previously done on the momentum equation, at least insofar as
    //the magnetic field is concerned.
    if(kinematic)
    {
        fillTimeField(u->vec, KINEMATIC);
    }

    //This is really the induction term, not advection, though it does contain
    //advection effects within it.
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

    //The 1 means we store the result in the force vectors.
    decomposeSolenoidal(B->sol, rhs, 1);
    debug("Magnetic forces done\n");
}

/*
 * This routine is in charge of the Temperature equation.  Virtually all
 * of the terms can be enabled or disabled by parameters read in through the
 * configuration file.  This equation has the form:
 * 
 * dT/dt = div(uT) + w hat z + del^2 T 
 */
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
        //advect the background profile (as long as it is enabled)
        if(tempBackground)
        {
            plusEq(forces, u->vec->z->spectral);
        }

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
    
    //Apply hyper diffusion to the boundaries
    if(sanitize)
    {
        killBoundaries(T->spatial, forces, 1, 100);
    }
}

/*
 * These are trash vectors that we will use for intermediate calculations.
 */
void initPhysics()
{
    temp1 = newVector(SPEC | SPAT);
    temp2 = newVector(SPEC | SPAT);
    rhs = newVector(SPEC);
}

/*
 * Cleanup after the init method.
 */
void finalizePhysics()
{
    deleteVector(&temp1);
    deleteVector(&temp2);
    deleteVector(&rhs);
}

/*
 * This is an entry point to the time integration.  When things are fully
 * running we do a third level Adams-Bashforth scheme which requires knowing the
 * past three forcing evaluations.  Since these are not available intitially,
 * the first few steps are lower order while we ramp up.
 */
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

/*
 * Horribly basic explicit euler step.  We just add dt * force to all of our 
 * variables.
 * 
 * For divergence free variables, we do time integration on the poloidal and
 * toroidal scalars, rather than on the vector itself.  This makes things
 * slightly verbose, as we then have to manually track the horizontal means
 * as well.
 */
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

/*
 * See AB3Step for more comments.
 */
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

/*
 * This is our main integration routine.  AB methods rely on knowing the past
 * state of the system for multiple times in the past.  In essence, this method
 * interpolates a curve between these last known points, and then performs an
 * exact integration over this curve to proceed from the last time step to the
 * next.  The fact that we have variable time steps complicates the derived 
 * coefficients, as can be seen by c0, c1 and c2.
 */
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

//This function is only designed to work on 2D y-invariant simulations
//It is also highly experimental and not fully functional.
/*
 * This routine calculates the center of mass of B2, in order to keep 
 * concentrations of magnetic field centered in the domain.
 */
displacement displacementByCenter()
{
    int i;
    displacement ret;
    complex PRECISION dkz;
    complex PRECISION dkx;
    
    complex PRECISION * data = B->vec->y->spectral;
    PRECISION mean = 0;
    PRECISION weightedMeanX = 0;
    PRECISION weightedMeanZ = 0;
    
    //Only one processor has the data to find the z center of mass, as it only
    //uses the mean values of x
    if(my_kx->min == 0)
    {
        mean = __real__(data[0]);
		weightedMeanZ = zmx/2.*mean;
        for(i = 1; i < ndkz; i++)
        {
            dkz = dzFactor(i);
            weightedMeanZ += __real__(data[i] / dkz);
        }
    }
    
    //Everyone shares in doing part of the x center mass
    if(vrank == 0)
	{
		weightedMeanX = xmx/2.*mean;
    	for(i = 1; i < my_kx->width; i++)
    	{
    	    dkx = dxFactor(i);
    	    weightedMeanX += 2*__real__(data[ndkz*i] / dkx);
   		}
	}
	else
	{
    	for(i = 0; i < my_kx->width; i++)
    	{
    	    dkx = dxFactor(i);
    	    weightedMeanX += 2*__real__(data[ndkz*i] / dkx);
   		}
    }
    //Root need all the info for the scattered sum
    if(vrank == 0)
    {
        MPI_Reduce(MPI_IN_PLACE, &weightedMeanX, 1, MPI_PRECISION, MPI_SUM, 0, vcomm);
        ret.dx = (weightedMeanX/mean) - xmx/2.;
        ret.dy = 0;
        ret.dz = (weightedMeanZ/mean) - zmx/2.;
    }   
    else
        MPI_Reduce(&weightedMeanX, 0, 1, MPI_PRECISION, MPI_SUM, 0, vcomm);
        
    MPI_Bcast(&ret, sizeof(displacement), MPI_BYTE, 0, vcomm);
    
    return ret;
}
