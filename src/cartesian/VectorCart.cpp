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

#include "cartesian/FieldCart.h"
#include "cartesian/VectorCart.h"
#include "cartesian/SolenoidCart.h"
#include "cartesian/TensorCart.h"
#include "VariableFactory.h"

using namespace std;

VectorCart::VectorCart()
{
    
}

shared_ptr<VectorCart> VectorCart::getShared()
{
    return static_pointer_cast<VectorCart>(shared_from_this());
}

shared_ptr<Field> VectorCart::divergence()
{
    shared_ptr<Field> ret = VariableFactory::createField();
    Div * node = new Div(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> VectorCart::curl()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    Curl * node = new Curl(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> VectorCart::laplacian()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    Laplacian * node = new Laplacian(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Tensor> VectorCart::gradient()
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    Grad * node = new Grad(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Solenoid> VectorCart::decompose()
{
    shared_ptr<Solenoid> ret = VariableFactory::createSolenoid();
    Decompose * node = new Decompose(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Solenoid> VectorCart::decomposeCurl()
{
    shared_ptr<Solenoid> ret = VariableFactory::createSolenoid();
    DecomposeCurl * node = new DecomposeCurl(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

VectorCart::Div::Div(std::shared_ptr<VectorCart> v)
        : VectorCart::OperatorTier2(v->op, div) {}

VectorCart::Grad::Grad(std::shared_ptr<VectorCart> v)
        : VectorCart::OperatorTier2(v->op, grad) {}

VectorCart::Curl::Curl(std::shared_ptr<VectorCart> v)
        : VectorCart::OperatorTier2(v->op, curl) {}

VectorCart::Laplacian::Laplacian(std::shared_ptr<VectorCart> v)
        : VectorCart::OperatorTier2(v->op, laplace) {}

VectorCart::Decompose::Decompose(std::shared_ptr<VectorCart> v)
        : VectorCart::OperatorTier1(v->op, decompose) {}

VectorCart::DecomposeCurl::DecomposeCurl(std::shared_ptr<VectorCart> v)
        : VectorCart::OperatorTier1(v->op, decomposeCurl) {}