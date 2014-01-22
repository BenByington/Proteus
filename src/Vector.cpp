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


#include "Vector.h"
#include "Scalar.h"
#include "Field2.h"
#include "Scalar.h"
#include "VariableFactory.h"

using namespace std;

Vector::Vector()
{
}

shared_ptr<Vector> Vector::multiply(shared_ptr<Scalar> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier3<Scalar,Vector> * node = new OperatorTier3<Scalar,Vector>(fact, shared_from_this(),OperatorTier3<Scalar,Vector>::mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier3<Scalar,Vector> * node = new OperatorTier3<Scalar,Vector>(fact, shared_from_this(),OperatorTier3<Scalar,Vector>::divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::multiply(shared_ptr<Field> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier3<Field,Vector> * node = new OperatorTier3<Field,Vector>(fact, shared_from_this(),OperatorTier3<Field,Vector>::mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::divide(shared_ptr<Field> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier3<Field,Vector> * node = new OperatorTier3<Field,Vector>(fact, shared_from_this(),OperatorTier3<Field,Vector>::divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::add(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier5<Vector,Vector> * node = new OperatorTier5<Vector,Vector>(shared_from_this(),r,OperatorTier5<Vector,Vector>::add);

    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Vector> Vector::subtract(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier5<Vector,Vector> * node = new OperatorTier5<Vector,Vector>(shared_from_this(),r,OperatorTier5<Vector,Vector>::sub);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Vector> Vector::cross(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier4<Vector,Vector> * node = new OperatorTier4<Vector,Vector>(shared_from_this(),r,OperatorTier4<Vector,Vector>::cross);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Vector::dot(shared_ptr<Vector> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier3<Vector,Vector> * node = new OperatorTier3<Vector,Vector>(shared_from_this(),r,OperatorTier3<Vector,Vector>::dot);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

shared_ptr<Tensor> Vector::outter(shared_ptr<Vector> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    OperatorTier4<Vector,Vector> * node = new OperatorTier4<Vector,Vector>(shared_from_this(),r,OperatorTier4<Vector,Vector>::outer);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}
