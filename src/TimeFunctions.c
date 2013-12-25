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

#include "TimeFunctions.h"
#include "Environment.h"
#include "Log.h"
#include "Communication.h"

#include <math.h>

void fillTimeField(p_vector vec, int func)
{
    int i,j,k;
    int index = 0;
    PRECISION x, y, z;
    if(func == MOMENTUM)
    {
        debug("Doing time dependant forcing of momentum equation\n");
        for(i = 0; i < my_z->width; i++)
        {
            z = getz(i);
            for(j = 0; j < my_x->width; j++)
            {
                x = getx(j);
                for(k = 0; k < ny; k++)
                {
                    y = gety(k);

                    vec->x->spatial[index] = mSourceX(x, y, z, elapsedTime);
                    vec->y->spatial[index] = mSourceY(x, y, z, elapsedTime);
                    vec->z->spatial[index] = mSourceZ(x, y, z, elapsedTime);

                    index++;
                }
            }
        }
    }
    else if(func == MAGNETIC)
    {
        debug("Doing time dependant forcing of induction equation\n");
        for(i = 0; i < my_z->width; i++)
        {
            z = getz(i);
            for(j = 0; j < my_x->width; j++)
            {
                x = getx(j);
                for(k = 0; k < ny; k++)
                {
                    y = gety(k);

                    vec->x->spatial[index] = bSourceX(x, y, z, elapsedTime);
                    vec->y->spatial[index] = bSourceY(x, y, z, elapsedTime);
                    vec->z->spatial[index] = bSourceZ(x, y, z, elapsedTime);

                    index++;
                }
            }
        }
    }
    else if(func == KINEMATIC)
    {
        debug("Doing time dependant velocity for kinematic problem\n");
        for(i = 0; i < my_z->width; i++)
        {
            z = getz(i);
            for(j = 0; j < my_x->width; j++)
            {
                x = getx(j);
                for(k = 0; k < ny; k++)
                {
                    y = gety(k);

                    vec->x->spatial[index] = kinematicX(x, y, z, elapsedTime);
                    vec->y->spatial[index] = kinematicY(x, y, z, elapsedTime);
                    vec->z->spatial[index] = kinematicZ(x, y, z, elapsedTime);

                    index++;
                }
            }
        }
    }
    else
    {
        error("ABORTING: Invalid function identifier: %d\n", func);
        return;
    }

    fftForward(vec->x);
    fftForward(vec->y);
    fftForward(vec->z);
}

extern PRECISION getx(int i)
{
    return 2 * PI * (PRECISION) (i + my_x->min) / (PRECISION) nx;
}

extern PRECISION gety(int i)
{
    return 2 * PI * (PRECISION) i  / (PRECISION) ny;
}

extern PRECISION getz(int i)
{
    return 2 * PI * (PRECISION) (i + my_z->min) / (PRECISION) nz;
}

/*
 * Kinematic flow field for the Brummel et al modified ABC flow
 */
inline PRECISION kinematicX(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
    return sin(z + momEps*sin(momOmega*t)) + cos(y + momEps*sin(momOmega*t));
}

inline PRECISION kinematicY(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
  return sin(x + momEps*sin(momOmega*t)) + cos(z + momEps*sin(momOmega*t));
}

inline PRECISION kinematicZ(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{

  return sin(y + momEps*sin(momOmega*t)) + cos(x + momEps*sin(momOmega*t));
}

/**
 * Forcing function that should result in the modified ABC flow
 **/

inline PRECISION mSourceX(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
    PRECISION temp1 = 0;
    PRECISION temp2 = 0;
    
    temp1 = cos(z + momEps*sin(momOmega*t)) - sin(y + momEps*sin(momOmega*t));
    temp2 = sin(z + momEps*sin(momOmega*t)) + cos(y + momEps*sin(momOmega*t));
    
    return momEps*momOmega*cos(momOmega*t) * temp1 + temp2 * Pr;
}

inline PRECISION mSourceY(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
    PRECISION temp1 = 0;
    PRECISION temp2 = 0;

    temp1 = cos(x + momEps*sin(momOmega*t)) - sin(z + momEps*sin(momOmega*t));
    temp2 = sin(x + momEps*sin(momOmega*t)) + cos(z + momEps*sin(momOmega*t));

    return momEps*momOmega*cos(momOmega*t) * temp1 + temp2 * Pr;
}

inline PRECISION mSourceZ(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
    PRECISION temp1 = 0;
    PRECISION temp2 = 0;

    temp1 = cos(y + momEps*sin(momOmega*t)) - sin(x + momEps*sin(momOmega*t));
    temp2 = sin(y + momEps*sin(momOmega*t)) + cos(x + momEps*sin(momOmega*t));

    return momEps*momOmega*cos(momOmega*t) * temp1 + temp2 * Pr;
}


/**
 * !!!MAKE SURE YOU PROGRAM IN A SOLENOIDAL FIELD!!!
 * If you do not, it will be projected into one in a potentially
 * unpredictable manor.
 *
 * Force induction equation in a similar manner to how
 * the momentum equation is forced.
 * Apply forcing func to B0sin(kz)sin(wt)[sin(ky),0,0]
 *
 * yields w*B0*sin(kz)*[sin(ky),0,0]*cos(wt)
 * + 2*k*k*Pr*B0*sin(kz)*[sin(ky),0,0]*sin(wt)/Pm
 *
 * Must pick suitable values for B0, k, w
 * Use Pr=0.01, Pm=1.0
 * Try B0=0.2, k=1.0, w=1e-4
 * should add field slower than decay rate
 * and such that it wont decay too rapidly
 *
 *
 * !!! possibly temporary change.  Removing the time dependence for the forcing
 * Apply forcing func to B0sin(kz)[sin(ky),0,0]
 *   yields 2*k*k*Pr*B0*sin(kz)*[sin(ky),0,0]/Pm
 **/
inline PRECISION bSourceX(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
  PRECISION temp1 = 0;
  PRECISION temp2 = 0;

  //temp1 = magW*magB0*sin(magK*z)*cos(magW*t)*sin(magK*y);
  //temp2 = 2.*magK*magK*magB0*sin(magK*z)*sin(magW*t)*sin(magK*y);

  temp1 = 0;
  temp2 = 2.*magK*magK*magB0*sin(magK*z)*sin(magK*y);

  return temp1 + temp2 * Pr / Pm;
}

inline PRECISION bSourceY(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
    return 0;
}

inline PRECISION bSourceZ(PRECISION x, PRECISION y, PRECISION z, PRECISION t)
{
    return 0;
}
