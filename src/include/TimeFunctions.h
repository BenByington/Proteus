/* 
 * File:   TimeFunctions.h
 * Author: Ben
 *
 * Created on July 9, 2010, 9:06 AM
 */

/** Usage Notes::
 *
 *  The user is to define the mSource, bSource and kinematic functions in order
 *  to get time dependant forcing for the momentum and magnetic equations, or
 *  a time dependant velocity in the kinematic problem.  Whatever function is
 *  defined uses a false coordinate system, and needs to satisfy f(0) = f(1)
 *  in all three dimensions.  The code will automatically stretch this to
 *  fit whatever the actual box dimensions are.  This simplifies things if you
 *  just want the forcing function to have a specific shape in the box, and
 *  complicates things if for whatever reason you are trying to specify the
 *  derivatives.  Also, when defining the functions keep in mind that funcions
 *  too complicated (too long, or doing fancy things) may cause the complier to
 *  ignore the 'inline' derictive.  The performance hit will likely be small or
 *  negligable, but it is still something to consider.
 **/

#ifndef _TIMEFUNCTIONS_H
#define	_TIMEFUNCTIONS_H

#define MOMENTUM 0
#define MAGNETIC 1
#define KINEMATIC 2

#include "Field.h"

void fillTimeField(p_vector, int func);

inline double getx(int i);
inline double gety(int j);
inline double getz(int j);

inline double mSourceX(double x, double y, double z, double t);
inline double mSourceY(double x, double y, double z, double t);
inline double mSourceZ(double x, double y, double z, double t);

inline double bSourceX(double x, double y, double z, double t);
inline double bSourceY(double x, double y, double z, double t);
inline double bSourceZ(double x, double y, double z, double t);

inline double kinematicX(double x, double y, double z, double t);
inline double kinematicY(double x, double y, double z, double t);
inline double kinematicZ(double x, double y, double z, double t);

#endif	/* _TIMEFUNCTIONS_H */

