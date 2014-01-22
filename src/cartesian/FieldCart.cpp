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
#include "VariableFactory.h"
using namespace std;

FieldCart::FieldCart()
{
    
}

shared_ptr<FieldCart> FieldCart::getShared()
{
    return static_pointer_cast<FieldCart>(shared_from_this());
}

shared_ptr<Field> FieldCart::laplacian()
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier2<Field> * node = new OperatorTier2<Field>(shared_from_this(),OperatorTier2<Field>::laplace);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

shared_ptr<Vector> FieldCart::gradient()
{
    shared_ptr<Vector> ret = VariableFactory::createVector();
    OperatorTier2<Field> * node = new OperatorTier2<Field>(shared_from_this(),OperatorTier2<Field>::grad);
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}
