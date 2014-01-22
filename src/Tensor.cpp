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


#include "Tensor.h"
#include "Vector.h"
#include "Scalar.h"
#include "VariableFactory.h"

using namespace std;

Tensor::Tensor()
{
    op = 0;
}

shared_ptr<Tensor> Tensor::multiply(shared_ptr<Scalar> fact)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    ScalarMul * node = new ScalarMul(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    ScalarDiv * node = new ScalarDiv(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::add(shared_ptr<Tensor> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    TensorAdd * node = new TensorAdd(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::subtract(shared_ptr<Tensor> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    TensorSub * node = new TensorSub(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Tensor::ScalarMul::ScalarMul(std::shared_ptr<Scalar> s, std::shared_ptr<Tensor> t)
        : OperatorTier3(t->op, s->op, mul) {}

Tensor::ScalarDiv::ScalarDiv(std::shared_ptr<Scalar> s, std::shared_ptr<Tensor> t)
        : OperatorTier3(t->op, s->op, divide) {}

Tensor::TensorAdd::TensorAdd(std::shared_ptr<Tensor> t1, std::shared_ptr<Tensor> t2)
        : OperatorTier5(t1->op, t2->op, add) {}

Tensor::TensorSub::TensorSub(std::shared_ptr<Tensor> t1, std::shared_ptr<Tensor> t2)
        : OperatorTier5(t1->op, t2->op, sub) {}