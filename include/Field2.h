/*
 * Copyright 2013-2014 Benjamin Byington
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

#ifndef FIELD2_H
#define FIELD2_H

#include "Scalar.h"
#include "Variable.h"

#include <memory>

class Vector;

class Field : public Variable, public std::enable_shared_from_this<Field>
{
protected:
    Field();
public:
    virtual ~Field(){}
    
    /*Apply a factor to our variable field*/
    std::shared_ptr<Field> multiply(std::shared_ptr<Scalar> fact);
    std::shared_ptr<Field> divide(std::shared_ptr<Scalar> mult);
    
    /*Arithmetic operations between variable fields*/
    std::shared_ptr<Field> add(std::shared_ptr<Field> r);
    std::shared_ptr<Field> subtract(std::shared_ptr<Field> r);
    std::shared_ptr<Field> multiply(std::shared_ptr<Field> r);
    std::shared_ptr<Field> divide(std::shared_ptr<Field> r);
    
    virtual std::shared_ptr<Vector> gradient() = 0;
    virtual std::shared_ptr<Field> laplacian() = 0;

private:
    class ScalarFactor : public GNode
    {
        friend class Field;
    public:
        ScalarFactor(std::shared_ptr<Scalar> s, std::shared_ptr<Field> v);
        virtual void execute();
        virtual std::string executeText();
        
    private:
        std::shared_ptr<Scalar> sParent;
        std::shared_ptr<Field> vParent;
        
        enum operations {mul, divide};
        
        operations op;
    };
    
    class FieldArithmetic : public GNode
    {
        friend class Field;
    public:
        FieldArithmetic(std::shared_ptr<Field> p1, std::shared_ptr<Field> p2);
        virtual void execute();
        virtual std::string executeText();
    private:
        std::shared_ptr<Field> p1;
        std::shared_ptr<Field> p2;
        
        enum operations {add, sub, mul, divide};
        operations op;
    };
};

#endif