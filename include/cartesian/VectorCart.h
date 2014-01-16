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

#ifndef VECTOR_CARTESIAN_H
#define VECTOR_CARTESIAN_H

#include "Vector.h"

class VectorCart : public Vector
{
    friend class FieldCart;
protected:
    VectorCart();
        
public:
    virtual ~VectorCart(){}
    
    virtual Solenoid * decompose();
    virtual Solenoid * decomposeCurl();
    virtual Field * divergence();
    virtual Vector * curl();

private:
    class AgnosticDeriv : public GNode
    {
        friend class VectorCart;
    public:
        AgnosticDeriv(VectorCart * v);
        virtual void execute();
        
    private:
        VectorCart * vParent;
        
        enum operations {div, curl};
        operations op;
    };
    
    class SolenoidOp : public GNode
    {
        friend class VectorCart;
    public:
        SolenoidOp(VectorCart * v);
        virtual void execute();
        
    private:
        VectorCart * vParent;
        
        enum operations {decomp, decompCurl};
        operations op;
    };
    
};
#endif