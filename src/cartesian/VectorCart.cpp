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

shared_ptr<Field> VectorCart::divergence()
{
    shared_ptr<Field> ret = VariableFactory::createField();
    AgnosticDeriv * node = new AgnosticDeriv(getShared());
    node->op = node->div;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> VectorCart::curl()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    AgnosticDeriv * node = new AgnosticDeriv(getShared());
    node->op = node->curl;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> VectorCart::laplacian()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    AgnosticDeriv * node = new AgnosticDeriv(getShared());
    node->op = node->laplace;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Tensor> VectorCart::gradient()
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    AgnosticDeriv * node = new AgnosticDeriv(getShared());
    node->op = node->grad;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Solenoid> VectorCart::decompose()
{
    shared_ptr<Solenoid> ret = VariableFactory::createSolenoid();
    AgnosticDeriv * node = new AgnosticDeriv(getShared());
    node->op = node->decomp;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Solenoid> VectorCart::decomposeCurl()
{
    shared_ptr<Solenoid> ret = VariableFactory::createSolenoid();
    AgnosticDeriv * node = new AgnosticDeriv(getShared());
    node->op = node->decompCurl;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

VectorCart::AgnosticDeriv::AgnosticDeriv(shared_ptr<VectorCart> v)
{
    this->vParent = v;
}

void VectorCart::AgnosticDeriv::execute()
{
    
}

string VectorCart::AgnosticDeriv::getDependString()
{
    string opName;
    switch(op)
    {
    case curl:
        opName = "Curl";
        break;
    case div:
        opName = "Divergence";
        break;
    case grad:
        opName = "Gradient";
        break;
    case decomp:
        opName = "Decompose Vector";
        break;
    case decompCurl:
        opName = "Decompose Curled Vector";
        break;
    }
    
    string ret = opName + ": " + this->vParent->op->getID();
    
    return ret;
}

shared_ptr<VectorCart> VectorCart::getShared()
{
    return static_pointer_cast<VectorCart>(shared_from_this());
}