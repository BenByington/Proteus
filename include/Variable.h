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

#ifndef VARIABLE_H
#define VARIABLE_H

#include "GNode.h"

class Vector;
class Scalar;

class Variable
{
protected:
    Variable();
    virtual Variable * createVariable() = 0;
    
public:
    virtual ~Variable(){}
    
    /*Apply a factor to our variable field*/
    Variable * operator *(Scalar * fact);
    Variable * operator /(Scalar * mult);
    
    /*Arithmetic operations between variable fields*/
    Variable * operator +(Variable * r);
    Variable * operator -(Variable * r);
    Variable * operator *(Variable * r);
    Variable * operator /(Variable * r);
    
    virtual Vector * gradient() = 0;
    virtual Variable * laplacian() = 0;
    
    GNode * op;

private:
    class ScalarFactor : public GNode
    {
        friend class Variable;
    public:
        ScalarFactor(Scalar * s, Variable * v);
        virtual void execute();
        
    private:
        Scalar * sParent;
        Variable * vParent;
        
        enum operations {mul, divide};
        operations op;
    };
    
    class VariableArithmetic : public GNode
    {
        friend class Variable;
    public:
        VariableArithmetic(Variable * p1, Variable * p2);
        virtual void execute();
    private:
        Variable * p1;
        Variable * p2;
        
        enum operations {add, sub, mul, divide};
        operations op;
    };
};

#endif