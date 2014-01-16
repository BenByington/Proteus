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

#ifndef FIELD2_H
#define FIELD2_H

#include "Scalar.h"

class Vector;

class Field
{
protected:
    Field();
public:
    virtual ~Field(){}
    virtual Field * createField() = 0;
    virtual Vector * createVector() = 0;
    
    /*Apply a factor to our variable field*/
    Field * operator *(Scalar * fact);
    Field * operator /(Scalar * mult);
    
    /*Arithmetic operations between variable fields*/
    Field * operator +(Field * r);
    Field * operator -(Field * r);
    Field * operator *(Field * r);
    Field * operator /(Field * r);
    
    virtual Vector * gradient() = 0;
    virtual Field * laplacian() = 0;
    
    GNode * op;

private:
    class ScalarFactor : public GNode
    {
        friend class Field;
    public:
        ScalarFactor(Scalar * s, Field * v);
        virtual void execute();
        virtual std::string executeText();
        
    private:
        Scalar * sParent;
        Field * vParent;
        
        enum operations {mul, divide};
        
        operations op;
    };
    
    class FieldArithmetic : public GNode
    {
        friend class Field;
    public:
        FieldArithmetic(Field * p1, Field * p2);
        virtual void execute();
        virtual std::string executeText();
    private:
        Field * p1;
        Field * p2;
        
        enum operations {add, sub, mul, divide};
        operations op;
    };
};

#endif