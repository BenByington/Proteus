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
        debug("Doing time dependant forcing of momentum equation",0);
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
        debug("Doing time dependant forcing of induction equation",0);
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
        debug("Doing time dependant velocity for kinematic problem",0);
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
        error("ABORTING: Invalid function identifier: %d", func);
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

inline double mSourceX(double x, double y, double z, double t)
{
    return 0;
}

inline double mSourceY(double x, double y, double z, double t)
{
    return 0;
}

inline double mSourceZ(double x, double y, double z, double t)
{
    return 0;
}

/**
 * !!!MAKE SURE YOU PROGRAM IN A SOLENOIDAL FIELD!!!
 * If you don't, strange and unpredictable things may happen
 * when it is decomposed into P and T...
 **/
inline double bSourceX(double x, double y, double z, double t)
{
    return 0;
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
    return 0;
}

inline double kinematicY(double x, double y, double z, double t)
{
    return 0;
}

inline double kinematicZ(double x, double y, double z, double t)
{
    return 0;
}
