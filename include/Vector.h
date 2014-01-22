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

#ifndef VECTOR_H
#define VECTOR_H

#include "Variable.h"

#include <string>
#include <memory>

class Field;
class Solenoid;
class Tensor;
class Scalar;

class Vector : public Variable, public std::enable_shared_from_this<Vector>
{
protected:
    Vector();
    
public:
    virtual ~Vector(){}
    
    std::shared_ptr<Vector> multiply(std::shared_ptr<Scalar> fact);
    std::shared_ptr<Vector> divide(std::shared_ptr<Scalar> mult);
    
    std::shared_ptr<Vector> multiply(std::shared_ptr<Field> f);
    std::shared_ptr<Vector> divide(std::shared_ptr<Field> f);
    
    /*Arithmetic operations between vectors fields*/
    std::shared_ptr<Vector> add(std::shared_ptr<Vector> r);
    std::shared_ptr<Vector> subtract(std::shared_ptr<Vector> r);
    
    std::shared_ptr<Field> dot(std::shared_ptr<Vector> r);
    std::shared_ptr<Vector> cross(std::shared_ptr<Vector> r);
    std::shared_ptr<Tensor> outter(std::shared_ptr<Vector> r);
    
    virtual std::shared_ptr<Solenoid> decompose() = 0;
    virtual std::shared_ptr<Solenoid> decomposeCurl() = 0;
    virtual std::shared_ptr<Field> divergence() = 0;
    virtual std::shared_ptr<Vector> curl() = 0;
    virtual std::shared_ptr<Tensor> gradient() = 0;
    virtual std::shared_ptr<Vector> laplacian() = 0;

};

#endif