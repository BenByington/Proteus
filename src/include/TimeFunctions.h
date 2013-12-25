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
#include "Precision.h"

void fillTimeField(p_vector, int func);

inline PRECISION getx(int i);
inline PRECISION gety(int j);
inline PRECISION getz(int j);

inline PRECISION mSourceX(PRECISION x, PRECISION y, PRECISION z, PRECISION t);
inline PRECISION mSourceY(PRECISION x, PRECISION y, PRECISION z, PRECISION t);
inline PRECISION mSourceZ(PRECISION x, PRECISION y, PRECISION z, PRECISION t);

inline PRECISION bSourceX(PRECISION x, PRECISION y, PRECISION z, PRECISION t);
inline PRECISION bSourceY(PRECISION x, PRECISION y, PRECISION z, PRECISION t);
inline PRECISION bSourceZ(PRECISION x, PRECISION y, PRECISION z, PRECISION t);

inline PRECISION kinematicX(PRECISION x, PRECISION y, PRECISION z, PRECISION t);
inline PRECISION kinematicY(PRECISION x, PRECISION y, PRECISION z, PRECISION t);
inline PRECISION kinematicZ(PRECISION x, PRECISION y, PRECISION z, PRECISION t);

#endif	/* _TIMEFUNCTIONS_H */

