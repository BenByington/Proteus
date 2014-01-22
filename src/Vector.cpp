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
    ScalarMul * node = new ScalarMul(fact, shared_from_this());
//    node->setOp(node->mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    ScalarDiv * node = new ScalarDiv(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::multiply(shared_ptr<Field> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    FieldMul * node = new FieldMul(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::divide(shared_ptr<Field> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    FieldDiv * node = new FieldDiv(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::add(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    VectorAdd * node = new VectorAdd(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Vector> Vector::subtract(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    VectorSub * node = new VectorSub(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Vector> Vector::cross(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    VectorCross * node = new VectorCross(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Vector::dot(shared_ptr<Vector> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    VectorDot * node = new VectorDot(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

shared_ptr<Tensor> Vector::outter(shared_ptr<Vector> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    VectorOuter * node = new VectorOuter(shared_from_this(), r);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

Vector::ScalarMul::ScalarMul(std::shared_ptr<Scalar> s, std::shared_ptr<Vector> v)
       : OperatorTier3(s->op, v->op, mul) {}

Vector::ScalarDiv::ScalarDiv(std::shared_ptr<Scalar> s, std::shared_ptr<Vector> v)
       : OperatorTier3(s->op, v->op, divide) {}

Vector::FieldMul::FieldMul(std::shared_ptr<Field> f, std::shared_ptr<Vector> v)
       : OperatorTier3(f->op, v->op, mul) {}

Vector::FieldDiv::FieldDiv(std::shared_ptr<Field> f, std::shared_ptr<Vector> v)
       : OperatorTier3(f->op, v->op, divide) {}

Vector::VectorAdd::VectorAdd(std::shared_ptr<Vector> v1, std::shared_ptr<Vector> v2)
       : OperatorTier5(v1->op, v2->op, add) {}

Vector::VectorSub::VectorSub(std::shared_ptr<Vector> v1, std::shared_ptr<Vector> v2)
       : OperatorTier5(v1->op, v2->op, sub) {}

Vector::VectorCross::VectorCross(std::shared_ptr<Vector> v1, std::shared_ptr<Vector> v2)
       : OperatorTier4(v1->op, v2->op, cross) {}

Vector::VectorOuter::VectorOuter(std::shared_ptr<Vector> v1, std::shared_ptr<Vector> v2)
       : OperatorTier4(v1->op, v2->op, outer) {}

Vector::VectorDot::VectorDot(std::shared_ptr<Vector> v1, std::shared_ptr<Vector> v2)
       : OperatorTier3(v1->op, v2->op, dot) {}