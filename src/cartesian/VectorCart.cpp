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

#include "cartesian/FieldCart.h"
#include "cartesian/VectorCart.h"
#include "cartesian/SolenoidCart.h"

using namespace std;

VectorCart::VectorCart()
{
    
}

Field * VectorCart::divergence()
{
    Field * ret = createField();
    AgnosticDeriv * node = new AgnosticDeriv(this);
    node->op = node->div;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

Vector * VectorCart::curl()
{
    Vector * ret = createVector();
    AgnosticDeriv * node = new AgnosticDeriv(this);
    node->op = node->curl;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

Solenoid * VectorCart::decompose()
{
    Solenoid * ret = createSolenoid();
    SolenoidOp * node = new SolenoidOp(this);
    node->op = node->decomp;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

Solenoid * VectorCart::decomposeCurl()
{
    Solenoid * ret = createSolenoid();
    SolenoidOp * node = new SolenoidOp(this);
    node->op = node->decompCurl;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

VectorCart::AgnosticDeriv::AgnosticDeriv(VectorCart* v)
{
    this->vParent = v;
}

void VectorCart::AgnosticDeriv::execute()
{
    
}

string VectorCart::AgnosticDeriv::executeText()
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
    }
    
    string ret = getName() + " = "; 
    ret += opName + ": " + this->vParent->op->getName();
    
    return ret;
}

VectorCart::SolenoidOp::SolenoidOp(VectorCart* v)
{
    this->vParent = v;
}

void VectorCart::SolenoidOp::execute()
{
    
}

string VectorCart::SolenoidOp::executeText()
{
    string opName;
    switch(op)
    {
    case decomp:
        opName = "Decompose Vector";
        break;
    case decompCurl:
        opName = "Decompose Curled Vector";
        break;
    }
    
    string ret = getName() + " = "; 
    ret += opName + ": " + this->vParent->op->getName();
    
    return ret;
}