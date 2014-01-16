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


#include "Tensor.h"
#include "Vector.h"
#include "Scalar.h"
#include "VariableFactory.h"

using namespace std;

Tensor::Tensor()
{
    op = 0;
}

Tensor * Tensor::operator *(Scalar * fact)
{
    Tensor * ret = VariableFactory::createTensor();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Tensor * Tensor::operator /(Scalar * fact)
{
    Tensor * ret = VariableFactory::createTensor();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Tensor * Tensor::operator +(Tensor * r)
{
    Tensor * ret = VariableFactory::createTensor();
    TensorArithmetic * node = new TensorArithmetic(this, r);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Tensor * Tensor::operator -(Tensor * r)
{
    Tensor * ret = VariableFactory::createTensor();
    TensorArithmetic * node = new TensorArithmetic(this, r);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}


Tensor::TensorArithmetic::TensorArithmetic(Tensor* v1, Tensor* v2)
{
    this->p1 = v1;
    this->p2 = v2;
}

void Tensor::TensorArithmetic::execute()
{
    
}

string Tensor::TensorArithmetic::executeText()
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
    
    string ret = getName() + " = "; 
    ret += opName + ": " + this->p1->op->getName();
    ret += " " + this->p2->op->getName();
    
    return ret;
}

Tensor::ScalarFactor::ScalarFactor(Scalar* s, Tensor * v)
{
    sParent = s;
    vParent = v;
}

void Tensor::ScalarFactor::execute()
{
    
}

string Tensor::ScalarFactor::executeText()
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
    
    string ret = getName() + string(" = "); 
    ret += opName + string(": ") + this->sParent->op->getName();
    ret += string(" ") + this->vParent->op->getName();
    
    
    return ret;
}