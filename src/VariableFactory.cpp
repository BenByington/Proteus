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

#include "VariableFactory.h"
#include "cartesian/periodic/FieldPeriodic.h"
#include "cartesian/periodic/VectorPeriodic.h"
#include "cartesian/periodic/SolenoidPeriodic.h"
#include "cartesian/periodic/TensorPeriodic.h"

#include <memory>

using namespace std;

shared_ptr<Field> VariableFactory::createField()
{
    return shared_ptr<Field>(new FieldPeriodic());
}

shared_ptr<Vector> VariableFactory::createVector()
{
    return shared_ptr<Vector>(new VectorPeriodic());
}

shared_ptr<Solenoid> VariableFactory::createSolenoid()
{
    return shared_ptr<Solenoid>(new SolenoidPeriodic());
}

shared_ptr<Tensor> VariableFactory::createTensor()
{
    return shared_ptr<Tensor>(new TensorPeriodic());
}

shared_ptr<Scalar> VariableFactory::createScalar()
{
    return shared_ptr<Scalar>(new Scalar());
}


shared_ptr<Field> VariableFactory::declareField(string name)
{
    shared_ptr<Field> ret = shared_ptr<Field>(new FieldPeriodic());
    ret->op->setLabel(name);
    return ret;
}

shared_ptr<Vector> VariableFactory::declareVector(string name)
{
    shared_ptr<Vector> ret = shared_ptr<Vector>(new VectorPeriodic());
    ret->op->setLabel(name);
    return ret;
}

shared_ptr<Solenoid> VariableFactory::declareSolenoid(string name)
{
    shared_ptr<Solenoid> ret = shared_ptr<Solenoid>(new SolenoidPeriodic());
    ret->op->setLabel(name);
    return ret;
}

shared_ptr<Tensor> VariableFactory::declareTensor(string name)
{
    shared_ptr<Tensor> ret = shared_ptr<Tensor>(new TensorPeriodic());
    ret->op->setLabel(name);
    return ret;
}

shared_ptr<Scalar> VariableFactory::declareScalar(string name)
{
    shared_ptr<Scalar> ret = shared_ptr<Scalar>(new Scalar());
    ret->op->setLabel(name);
    return ret;
}