/*
 * Copyright 2013 Benjamin Byington
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
#include "ExecutionGraph.h"

#include <memory>

using namespace std;

TopDagsTest::TopDagsTest()
{
    
}

void TopDagsTest::setup()
{
    shared_ptr<Vector> B = VariableFactory::declareVector("B");
    shared_ptr<Vector> u = VariableFactory::declareVector("u");
    shared_ptr<Field> T = VariableFactory::declareField("T");

    shared_ptr<Vector> hatz = VariableFactory::declareVector("z_hat");
    
    shared_ptr<Scalar> Pr = VariableFactory::declareScalar("Pr");
    shared_ptr<Scalar> Pm = VariableFactory::declareScalar("Pm");
    shared_ptr<Scalar> Ra = VariableFactory::declareScalar("Ra");
    shared_ptr<Scalar> alpha = VariableFactory::declareScalar("alpha");
    shared_ptr<Scalar> Beta = VariableFactory::declareScalar("Beta");
    
    ExecutionGraph * g = new ExecutionGraph();
    g->registerHead(B->op);
    g->registerHead(u->op);
    g->registerHead(T->op);
    g->registerHead(Pr->op);
    g->registerHead(Pm->op);
    g->registerHead(Ra->op);
    g->registerHead(alpha->op);
    g->registerHead(Beta->op);
    g->registerHead(hatz->op);
    
    shared_ptr<Vector> db = ((u->cross(B))->curl())->add(B->laplacian()->multiply(Pr)->divide(Pm));
    shared_ptr<Vector> du = u->laplacian()->multiply(Pr);  //diffusion
    
    shared_ptr<Tensor> tens = B->outter(B)->multiply(alpha);
    tens = tens->subtract(u->outter(u));
    
    du = du->add(tens->divergence()); //advection and lorentz
    shared_ptr<Field> buoy = T->multiply(Pr)->multiply(Ra)->add(B->dot(B)->divide(Beta));
    du = du->add(hatz->multiply(buoy));
    
    shared_ptr<Field> dt = T->laplacian()->add(hatz->dot(u));
    dt = dt->add(u->multiply(T)->divergence());
    
    shared_ptr<Solenoid> bs = db->decompose();
    shared_ptr<Solenoid> us = du->decomposeCurl();
    
    g->execute();
}

void TopDagsTest::cleanup()
{
    
}

bool TopDagsTest::run()
{
    
}