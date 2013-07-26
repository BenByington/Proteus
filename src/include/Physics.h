/* 
 * File:   Physics.h
 * Author: Ben
 *
 * Created on April 15, 2010, 8:45 AM
 */

#ifndef _PHYSICS_H
#define	_PHYSICS_H

#include "Precision.h"
#include "Environment.h"

void iterate();

void initPhysics();
void finalizePhysics();

//functions that can be used to find displacement for field shifting
displacement displacementByCenter();

#endif	/* _PHYSICS_H */

