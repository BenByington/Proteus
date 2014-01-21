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
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::add(shared_ptr<Tensor> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    TensorArithmetic * node = new TensorArithmetic(shared_from_this(), r);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Tensor> Tensor::subtract(shared_ptr<Tensor> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    TensorArithmetic * node = new TensorArithmetic(shared_from_this(), r);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}


Tensor::TensorArithmetic::TensorArithmetic(shared_ptr<Tensor> v1, shared_ptr<Tensor> v2)
{
    this->p1 = v1;
    this->p2 = v2;
}

void Tensor::TensorArithmetic::execute()
{
    
}

string Tensor::TensorArithmetic::getDependString()
{
    string opName;
    switch(op)
    {
    case sub:
        opName = "Subtract";
        break;
    case add:
        opName = "Add";    
    }
    
    string ret = opName + ": " + this->p1->op->getID();
    ret += " " + this->p2->op->getID();
    
    return ret;
}

Tensor::ScalarFactor::ScalarFactor(shared_ptr<Scalar> s, shared_ptr<Tensor> v)
{
    sParent = s;
    vParent = v;
}

void Tensor::ScalarFactor::execute()
{
    
}

string Tensor::ScalarFactor::getDependString()
{
    string opName;
    switch(op)
    {
    case mul:
        opName = "Multiply";
        break;
    case divide:
        opName = "Divide";
        break;
    }
    
    string ret = opName + string(": ") + this->sParent->op->getID();
    ret += string(" ") + this->vParent->op->getID();
    
    
    return ret;
}