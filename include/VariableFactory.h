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
#ifndef VARIABLEFACTORY_H
#define	VARIABLEFACTORY_H

#include "Scalar.h"
#include "Field2.h"
#include "Vector.h"
#include "Solenoid.h"
#include "Tensor.h"

#include <memory>
#include <string>

class VariableFactory
{
private:
    VariableFactory();
  
public:
    static std::shared_ptr<Scalar> createScalar();
    static std::shared_ptr<Field> createField();
    static std::shared_ptr<Vector> createVector();
    static std::shared_ptr<Solenoid> createSolenoid();
    static std::shared_ptr<Tensor> createTensor();
    
    static std::shared_ptr<Scalar> declareScalar(std::string name);
    static std::shared_ptr<Field> declareField(std::string name);
    static std::shared_ptr<Vector> declareVector(std::string name);
    static std::shared_ptr<Solenoid> declareSolenoid(std::string name);
    static std::shared_ptr<Tensor> declareTensor(std::string name);
};

#endif	

