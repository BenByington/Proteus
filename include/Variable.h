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

#ifndef VARIABLE_H
#define VARIABLE_H

#include "GNode.h"

#include <memory>

class Vector;
class Scalar;

class Variable
{
protected:
    Variable(){op = new nop();}

    class MathNode : public GNode
    {
    public:
        virtual std::string executeText();
        virtual std::string getDependString() = 0;
        virtual void updateLabel() = 0;
        ~MathNode(){}
    protected:
        MathNode();
        
        std::string opName;
        std::string opSym;
    };
    
    class nop : public MathNode
    {
    public:
        nop(){};
        virtual ~nop(){}
        virtual void execute(){};
        virtual std::string getDependString(){return "";}
        virtual void updateLabel(){};
    };
    
    class UnaryOp : public MathNode
    {
    public:
        virtual std::string getDependString();
        virtual void updateLabel();
        virtual ~UnaryOp(){}
    protected:
        UnaryOp(GNode * a);
        GNode * r;
    private:
        UnaryOp();
    };
    
    class BinaryOp : public MathNode
    {
    public:
        virtual std::string getDependString();
        virtual void updateLabel();
        virtual ~BinaryOp(){}
    protected:
        BinaryOp(GNode * a, GNode * b);
        GNode * l;
        GNode * r;
    private:
        BinaryOp();
    };
    
    template <class T>
    class OperatorTier1 : public UnaryOp
    {
    public:
        enum operation {decompose, decomposeCurl, recompose};
        OperatorTier1(std::shared_ptr<T> left, operation op);
        virtual ~OperatorTier1(){};
        virtual void execute();
    private:
        OperatorTier1();
    };
    
    template <class T>
    class OperatorTier2 : public UnaryOp
    {
    public:
        enum operation {div, grad, curl, laplace};
        OperatorTier2(std::shared_ptr<T> left, operation op);
        virtual ~OperatorTier2(){};
        virtual void execute();
    private:
        OperatorTier2();
    };
    
    template <class T1, class T2>
    class OperatorTier3 : public BinaryOp
    {
    public:
        enum operation {mul, divide, dot};
        OperatorTier3(std::shared_ptr<T1> left, std::shared_ptr<T2> right, operation op);
        virtual ~OperatorTier3(){};
        virtual void execute();
    private:
        OperatorTier3();
    };
    
    template <class T1, class T2>
    class OperatorTier4 : public BinaryOp
    {
    public:
        enum operation {cross, outer};
        OperatorTier4(std::shared_ptr<T1> left, std::shared_ptr<T2> right, operation op);
        virtual ~OperatorTier4(){};
        virtual void execute();
    private:
        OperatorTier4();
    };
    
    template <class T1, class T2>
    class OperatorTier5 : public BinaryOp
    {
    public:
        enum operation {add, sub};
        OperatorTier5(std::shared_ptr<T1> left, std::shared_ptr<T2> right, operation op);
        virtual ~OperatorTier5(){};
        virtual void execute();
    private:
        OperatorTier5();
    };
    
public:
    virtual ~Variable(){}
    
    GNode * op;

};

#endif