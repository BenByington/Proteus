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

#include "cartesian/SolenoidCart.h"

SolenoidCart::SolenoidCart()
{
    
}

Vector * SolenoidCart::recompose()
{
    Vector * ret = createVector();
    VectorOp * node = new VectorOp(this);
    node->op = node->recompose;
    ret->op = node;
    
    node->addDependency(this->op);
    
    return ret;
}

SolenoidCart::VectorOp::VectorOp(SolenoidCart * s)
{
    this->vParent = s;
}

void SolenoidCart::VectorOp::execute()
{
    
}