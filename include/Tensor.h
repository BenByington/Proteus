/*
 * Copywrite 2013-2014 Benjamin Byington
 *
 * This file is part of the Proteus software package
 * 
 * Proteus is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Proteus is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with Proteus.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef TENSOR_H
#define	TENSOR_H

#include "Variable.h"

class Vector;
class Scalar;

class Tensor : public Variable
{
protected:
    Tensor();
    
public:
    virtual ~Tensor(){}
    
    Tensor * operator *(Scalar * fact);
    Tensor * operator /(Scalar * mult);
    
    /*Arithmetic operations between vectors fields*/
    Tensor * operator +(Tensor * r);
    Tensor * operator -(Tensor * r);
    
    virtual Vector * divergence() = 0;
    
private:
    class TensorArithmetic : public GNode
    {
        friend class Tensor;
    public:
        TensorArithmetic(Tensor * v1, Tensor * v2);
        virtual void execute();
        virtual std::string executeText();
        
    private:
        Tensor * p1;
        Tensor * p2;
        
        enum operations {add, sub};
        operations op;
    };
    
    class ScalarFactor : public GNode
    {
        friend class Tensor;
    public:
        ScalarFactor(Scalar * s, Tensor * v);
        virtual void execute();
        virtual std::string executeText();
        
    private:
        Scalar * sParent;
        Tensor * vParent;
        
        enum operations {mul, divide};        
        operations op;
    };
};

#endif	/* TENSOR_H */

