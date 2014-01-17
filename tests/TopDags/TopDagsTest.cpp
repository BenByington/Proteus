/*
 * Copywrite 2013 Benjamin Byington
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

#include "TopDags/TopDagsTest.h"

#include "VariableFactory.h"
#include "Field.h"
#include "Vector.h"
#include "Field.h"
#include "Solenoid.h"
#include "Scalar.h"

TopDagsTest::TopDagsTest()
{
    
}

void TopDagsTest::setup()
{
    Vector * B = VariableFactory::createVector();
    Vector * u = VariableFactory::createVector();
    Field * T = VariableFactory::createField();
    
    Scalar * Rm = new Scalar();
    Scalar * Re = new Scalar();
    Scalar * Beta = new Scalar();
    Scalar * Buoy = new Scalar();
    Scalar * kappa = new Scalar();
    
    Vector * db = ((u->cross(B))->curl())->operator +(B->laplacian()->operator *(Rm));
    Vector * du = u->laplacian()->operator *(Re);  //diffusion
    
    Tensor * tens = B->outter(B);
    tens = tens->operator -(u->outter(u));
    
    du = du->operator +(tens->divergence()); //advection and lorentz
    
    
    Vector * dT = T->laplacian() * kappa;
}

void TopDagsTest::cleanup()
{
    
}

bool TopDagsTest::run()
{
    
}