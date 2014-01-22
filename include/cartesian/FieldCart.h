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

#ifndef Field_CARTESIAN_H
#define Field_CARTESIAN_H

#include "Field2.h"

#include <string>

class FieldCart : public Field
{
protected:
    FieldCart();
    
public:
    virtual ~FieldCart(){}
    
    virtual std::shared_ptr<Vector> gradient();
    virtual std::shared_ptr<Field> laplacian();

private:
    std::shared_ptr<FieldCart> getShared();
    
    class Grad : public OperatorTier2
    {
    public:
        Grad(std::shared_ptr<FieldCart> v);
        virtual void execute() {}
    };
    
    class Laplacian : public OperatorTier2
    {
    public:
        Laplacian(std::shared_ptr<FieldCart> v);
        virtual void execute() {}
    };
};

#endif