/* 
 * File:   State.h
 * Author: Ben
 *
 * Created on March 11, 2010, 2:38 PM
 */

#include "Field.h"
#include "Precision.h"

#ifndef _STATE_H
#define	_STATE_H

void initState();
void finalizeState();

extern p_componentVar B;
extern p_componentVar u;
extern p_field T;

extern PRECISION maxVel[3];
extern p_field forceField;


#endif	/* _STATE_H */

