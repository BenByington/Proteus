/*
 * Copywrite 2013-2014 Benjamin Byington
 *
 * This file is part of the IMHD software package
 * 
 * IMHD is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
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

#ifndef VECTOR_H
#define VECTOR_H

#include "GNode.h"
#include "Variable.h"

class Vector : public Variable
{
public:
    virtual ~Vector(){}
    
    /*Arithmetic operations between vectors fields*/
    virtual Vector * operator +(const Vector * r) = 0;
    virtual Vector * operator -(const Vector * r) = 0;
    
    virtual Variable * dot(const Vector * r) = 0;
    virtual Vector * cross(const Vector * r) = 0;
    
    virtual Variable * divergence() = 0;
    virtual Vector * curl() = 0;
};

#endif