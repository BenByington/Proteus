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

#ifndef VECTOR_CARTESIAN_H
#define VECTOR_CARTESIAN_H

#include "Vector.h"

class VectorCart : public Vector
{
    //friend class FieldCart;
protected:
    VectorCart();
        
public:
    virtual ~VectorCart(){}
    
    virtual std::shared_ptr<Solenoid> decompose();
    virtual std::shared_ptr<Solenoid> decomposeCurl();
    virtual std::shared_ptr<Field> divergence();
    virtual std::shared_ptr<Vector> curl();
    virtual std::shared_ptr<Tensor> gradient();
    virtual std::shared_ptr<Vector> laplacian();

private:
    std::shared_ptr<VectorCart> getShared();
    
    class Div : public Vector::OperatorTier2
    {
    public:
        Div(std::shared_ptr<VectorCart> v);
        virtual void execute() {}
    };
    
    class Grad : public Vector::OperatorTier2
    {
    public:
        Grad(std::shared_ptr<VectorCart> v);
        virtual void execute() {}
    };
    
    class Curl : public Vector::OperatorTier2
    {
    public:
        Curl(std::shared_ptr<VectorCart> v);
        virtual void execute() {}
    };
    
    class Laplacian : public Vector::OperatorTier2
    {
    public:
        Laplacian(std::shared_ptr<VectorCart> v);
        virtual void execute() {}
    };
    
    class Decompose : public Vector::OperatorTier1
    {
    public:
        Decompose(std::shared_ptr<VectorCart> v);
        virtual void execute() {}
    };
    
    class DecomposeCurl : public Vector::OperatorTier1
    {
    public:
        DecomposeCurl(std::shared_ptr<VectorCart> v);
        virtual void execute() {}
    };
    
};
#endif