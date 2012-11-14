#include "TimeFunctions.h"
#include "Environment.h"
#include "Log.h"
#include "Communication.h"

#include <math.h>

void fillTimeField(p_vector vec, int func)
{
    int i,j,k;
    int index = 0;
    double x, y, z;
    if(func == MOMENTUM)
    {
        debug("Doing time dependant forcing of momentum equation\n",0);
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
        debug("Doing time dependant forcing of induction equation\n",0);
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
        debug("Doing time dependant velocity for kinematic problem\n",0);
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

extern double getx(int i)
{
    return (double) (i + my_x->min) / (double) nx;
}

extern double gety(int i)
{
    return (double) i  / (double) ny;
}

extern double getz(int i)
{
    return (double) (i + my_z->min) / (double) nz;
}

/**
 * Forcing function, F (see Eq. 5 of Brummell et al. 1998)
 * applied to modified ABC flow
 * use epsilon = 1.0
 * Re = 100.0
 * Omega = 2.5 (vary this to change nonlinear dynamo properties)
 **/

inline double mSourceX(double x, double y, double z, double t)
{
  return 1.0*cos(1.0*t)*(cos(2.0*PI*z+sin(1.0*t))-sin(2.0*PI*y+sin(1.0*t)))+(sin(2.0*PI*z+sin(1.0*t))+cos(2.0*PI*y+sin(1.0*t)))/100.0;
}

inline double mSourceY(double x, double y, double z, double t)
{
  return 1.0*cos(1.0*t)*(cos(2.0*PI*x+sin(1.0*t))-sin(2.0*PI*z+sin(1.0*t)))+(sin(2.0*PI*x+sin(1.0*t))+cos(2.0*PI*z+sin(1.0*t)))/100.0;
}

inline double mSourceZ(double x, double y, double z, double t)
{
  return 1.0*cos(1.0*t)*(cos(2.0*PI*y+sin(1.0*t))-sin(2.0*PI*x+sin(1.0*t)))+(sin(2.0*PI*y+sin(1.0*t))+cos(2.0*PI*x+sin(1.0*t)))/100.0;
}


/**
 * !!!MAKE SURE YOU PROGRAM IN A SOLENOIDAL FIELD!!!
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
 **/
inline double bSourceX(double x, double y, double z, double t)
{
  static double k = 1.0;
  static double w = 1e-2;
  static double b0 = 0.1;

  return w*b0*sin(2*PI*k*z)*sin(2*PI*k*y)*cos(w*t)+2.0*k*k*0.01*b0*sin(2*PI*k*z)*sin(2*PI*k*y)*sin(w*t)/1.0;
}

inline double bSourceY(double x, double y, double z, double t)
{
    return 0;
}

inline double bSourceZ(double x, double y, double z, double t)
{
    return 0;
}

inline double kinematicX(double x, double y, double z, double t)
{
    return sin(2.0*PI*z+1.0*sin(1.0*t)) + cos(2.0*PI*y+1.0*sin(1.0*t));
}

inline double kinematicY(double x, double y, double z, double t)
{
  return sin(2.0*PI*x+1.0*sin(1.0*t)) + cos(2.0*PI*z+1.0*sin(1.0*t));
}

inline double kinematicZ(double x, double y, double z, double t)
{

  return sin(2*PI*y+1.0*sin(1.0*t)) + cos(2.0*PI*x+1.0*sin(1.0*t));
}
