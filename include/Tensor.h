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

#ifndef TENSOR_H
#define	TENSOR_H

#include "Variable.h"

#include <memory>

class Vector;
class Scalar;

class Tensor : public Variable, public std::enable_shared_from_this<Tensor>
{
protected:
    Tensor();
    
public:
    virtual ~Tensor(){}
    
    std::shared_ptr<Tensor> multiply(std::shared_ptr<Scalar> fact);
    std::shared_ptr<Tensor> divide(std::shared_ptr<Scalar> mult);
    
    /*Arithmetic operations between vectors fields*/
    std::shared_ptr<Tensor> add(std::shared_ptr<Tensor> r);
    std::shared_ptr<Tensor> subtract(std::shared_ptr<Tensor> r);
    
    virtual std::shared_ptr<Vector> divergence() = 0;
    
private:
    class ScalarMul : public OperatorTier3
    {
    public:
        ScalarMul(std::shared_ptr<Scalar> s, std::shared_ptr<Tensor> t);
        virtual void execute() {}
    };
    
    class ScalarDiv : public OperatorTier3
    {
    public:
        ScalarDiv(std::shared_ptr<Scalar> s, std::shared_ptr<Tensor> t);
        virtual void execute() {}
    };

    class TensorAdd : public OperatorTier5
    {
    public:
        TensorAdd(std::shared_ptr<Tensor> t1, std::shared_ptr<Tensor> t2);
        virtual void execute() {}
    };
    
    class TensorSub : public OperatorTier5
    {
    public:
        TensorSub(std::shared_ptr<Tensor> t1, std::shared_ptr<Tensor> t2);
        virtual void execute() {}
    };
    
};

#endif	/* TENSOR_H */

