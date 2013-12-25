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

//field to be used for boundary hyper diffusion term
extern p_field hyper;
extern p_field hyperWork;

extern PRECISION maxVel[3];
extern p_field forceField;
extern p_field magForceField;


#endif	/* _STATE_H */

