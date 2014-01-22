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
    OperatorTier2<Vector> * node = new OperatorTier2<Vector>(shared_from_this(),OperatorTier2<Vector>::div);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> VectorCart::curl()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier2<Vector> * node = new OperatorTier2<Vector>(shared_from_this(),OperatorTier2<Vector>::curl);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> VectorCart::laplacian()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier2<Vector> * node = new OperatorTier2<Vector>(shared_from_this(),OperatorTier2<Vector>::laplace);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Tensor> VectorCart::gradient()
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    OperatorTier2<Vector> * node = new OperatorTier2<Vector>(shared_from_this(),OperatorTier2<Vector>::grad);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Solenoid> VectorCart::decompose()
{
    shared_ptr<Solenoid> ret = VariableFactory::createSolenoid();
    OperatorTier1<Vector> * node = new OperatorTier1<Vector>(shared_from_this(),OperatorTier1<Vector>::decompose);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Solenoid> VectorCart::decomposeCurl()
{
    shared_ptr<Solenoid> ret = VariableFactory::createSolenoid();
    OperatorTier1<Vector> * node = new OperatorTier1<Vector>(shared_from_this(),OperatorTier1<Vector>::decomposeCurl);ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}
