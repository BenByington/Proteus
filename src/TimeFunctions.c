#include "TimeFunctions.h"
#include "Environment.h"
#include "Log.h"

#include <math.h>

void fillTimeField(p_field field, int func)
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

                    field->spatial[index] = mSource(x, y, z, elapsedTime);

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

                    field->spatial[index] = bSource(x, y, z, elapsedTime);

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

                    field->spatial[index] = kinematic(x, y, z, elapsedTime);

                    index++;
                }
            }
        }
    }
    else
    {
        error("ABORTING: Invalid function identifier: %d", func)
    }
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

inline double mSource(double x, double y, double z, double t)
{
    return 0;
}

inline double bSource(double x, double y, double z, double t)
{
    return 0;
}

inline double kinematic(double x, double y, double z, double t)
{
    return 0;
}
