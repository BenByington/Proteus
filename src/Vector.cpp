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
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->setOp(node->mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->setOp(node->divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::multiply(shared_ptr<Field> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    FieldFactor * node = new FieldFactor(fact, shared_from_this());
    node->setOp(node->mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::divide(shared_ptr<Field> fact)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    FieldFactor * node = new FieldFactor(fact, shared_from_this());
    node->setOp(node->divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Vector> Vector::add(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    VectorArithmetic * node = new VectorArithmetic(shared_from_this(), r);
    node->setOp(node->add);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Vector> Vector::subtract(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    VectorArithmetic * node = new VectorArithmetic(shared_from_this(), r);
    node->setOp(node->sub);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Vector> Vector::cross(shared_ptr<Vector> r)
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    VectorArithmetic * node = new VectorArithmetic(shared_from_this(), r);
    node->setOp(node->cross);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Vector::dot(shared_ptr<Vector> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    VectorArithmetic * node = new VectorArithmetic(shared_from_this(), r);
    node->setOp(node->dot);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

shared_ptr<Tensor> Vector::outter(shared_ptr<Vector> r)
{
    shared_ptr<Tensor> ret = VariableFactory::createTensor();
    VectorArithmetic * node = new VectorArithmetic(shared_from_this(), r);
    node->setOp(node->outter);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

Vector::VectorArithmetic::VectorArithmetic(shared_ptr<Vector> v1, shared_ptr<Vector> v2)
{
    this->p1 = v1;
    this->p2 = v2;
}

void Vector::VectorArithmetic::execute()
{
    
}

string Vector::VectorArithmetic::getDependString()
{
    string opName;
    switch(op)
    {
    case cross:
        opName = "Cross";
        break;
    case dot:
        opName = "Dot";
        break;
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

void Vector::VectorArithmetic::setOp(operations o)
{
    this->op = o;
    
    string opName;
    switch(op)
    {
    case cross:
        opName = string("\u2A2f");
        break;
    case dot:
        opName = string("\u00B7");
        break;
    case sub:
        opName = " - ";
        break;
    case add:
        opName = " + ";    
    }
    
    label = "(" + this->p1->op->getLabel() + opName + this->p2->op->getLabel() + ")";
}

Vector::ScalarFactor::ScalarFactor(shared_ptr<Scalar> s, shared_ptr<Vector> v)
{
    sParent = s;
    vParent = v;
}

void Vector::ScalarFactor::execute()
{
    
}

string Vector::ScalarFactor::getDependString()
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

void Vector::ScalarFactor::setOp(operations o)
{
    this->op = o;
    
    string opName;
    switch(op)
    {
    case mul:
        opName = " * ";
        break;
    case divide:
        opName = " / ";
        break;
    }
    
    label = "(" + this->vParent->op->getLabel() + opName + this->sParent->op->getLabel() + ")";
}

Vector::FieldFactor::FieldFactor(shared_ptr<Field> s, shared_ptr<Vector> v)
{
    fParent = s;
    vParent = v;
}

void Vector::FieldFactor::execute()
{
    
}

string Vector::FieldFactor::getDependString()
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
    
    string ret = opName + string(": ") + this->fParent->op->getID();
    ret += string(" ") + this->vParent->op->getID();
    
    
    return ret;
}

void Vector::FieldFactor::setOp(operations o)
{
    this->op = o;
    
    string opName;
    switch(op)
    {
    case mul:
        opName = " * ";
        break;
    case divide:
        opName = " / ";
        break;
    }
    
    label = "(" + this->vParent->op->getLabel() + opName + this->fParent->op->getLabel() + ")";
}