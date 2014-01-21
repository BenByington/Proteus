
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

#include "cartesian/TensorCart.h"
#include "cartesian/VectorCart.h"
#include "VariableFactory.h"

using namespace std;

TensorCart::TensorCart()
{
    
}

shared_ptr<Vector> TensorCart::divergence()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    AgnosticDeriv * node = new AgnosticDeriv(this);
    node->op = node->div;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

TensorCart::AgnosticDeriv::AgnosticDeriv(TensorCart* v)
{
    this->vParent = v;
}

void TensorCart::AgnosticDeriv::execute()
{
    
}

string TensorCart::AgnosticDeriv::getDependString()
{
    string opName = "Divergence";
    
    string ret = opName + ": " + this->vParent->op->getID();
    
    return ret;
}
