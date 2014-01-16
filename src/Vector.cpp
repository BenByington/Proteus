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


#include "Vector.h"
#include "Scalar.h"
#include "Field2.h"
#include "Scalar.h"
#include "VariableFactory.h"

using namespace std;

Vector::Vector()
{
    op = 0;
}

Vector * Vector::operator *(Scalar * fact)
{
    Vector * ret = VariableFactory::createVector();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Vector * Vector::operator /(Scalar * fact)
{
    Vector * ret = VariableFactory::createVector();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Vector * Vector::operator +(Vector * r)
{
    Vector * ret = VariableFactory::createVector();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Vector * Vector::operator -(Vector * r)
{
    Vector * ret = VariableFactory::createVector();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Vector * Vector::cross(Vector * r)
{
    Vector * ret = VariableFactory::createVector();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->cross;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Vector::dot(Vector* r)
{
    Field * ret = VariableFactory::createField();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->dot;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

Vector::VectorArithmetic::VectorArithmetic(Vector* v1, Vector* v2)
{
    this->p1 = v1;
    this->p2 = v2;
}

void Vector::VectorArithmetic::execute()
{
    
}

string Vector::VectorArithmetic::executeText()
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
    
    string ret = getName() + " = "; 
    ret += opName + ": " + this->p1->op->getName();
    ret += " " + this->p2->op->getName();
    
    return ret;
}

Vector::ScalarFactor::ScalarFactor(Scalar* s, Vector * v)
{
    sParent = s;
    vParent = v;
}

void Vector::ScalarFactor::execute()
{
    
}

string Vector::ScalarFactor::executeText()
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