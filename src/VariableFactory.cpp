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
