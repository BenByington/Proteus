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
    OperatorTier3<Scalar,Tensor> * node = new OperatorTier3<Scalar,Tensor>(fact,shared_from_this(),OperatorTier3<Scalar,Tensor>::mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    OperatorTier3<Scalar,Tensor> * node = new OperatorTier3<Scalar,Tensor>(fact,shared_from_this(),OperatorTier3<Scalar,Tensor>::divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::add(shared_ptr<Tensor> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    OperatorTier5<Tensor,Tensor> * node = new OperatorTier5<Tensor,Tensor>(shared_from_this(),r,OperatorTier5<Tensor,Tensor>::add);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::subtract(shared_ptr<Tensor> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    OperatorTier5<Tensor,Tensor> * node = new OperatorTier5<Tensor,Tensor>(shared_from_this(),r,OperatorTier5<Tensor,Tensor>::sub);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}
