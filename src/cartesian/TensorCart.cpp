
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

shared_ptr<TensorCart> TensorCart::getShared()
{
    return static_pointer_cast<TensorCart>(shared_from_this());
}

shared_ptr<Vector> TensorCart::divergence()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    Div * node = new Div(getShared());
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

TensorCart::Div::Div(std::shared_ptr<TensorCart> v)
        : OperatorTier2(v->op, div) {}