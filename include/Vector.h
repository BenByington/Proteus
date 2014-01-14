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
protected:
    Vector();
    virtual Vector * createVector() = 0;
    
public:
    virtual ~Vector(){}
    
    /*Arithmetic operations between vectors fields*/
    Vector * operator +(Vector * r);
    Vector * operator -(Vector * r);
    
    Variable * dot(Vector * r);
    Vector * cross(Vector * r);
    
    virtual Variable * divergence() = 0;
    virtual Vector * curl() = 0;

private:
    class VectorArithmetic : public GNode
    {
        friend class Vector;
    public:
        VectorArithmetic(Vector * v1, Vector * v2);
        virtual void execute();
        
    private:
        Vector * p1;
        Vector * p2;
        
        enum operations {add, sub, dot, cross};
        operations op;
    };
    
    
};

#endif