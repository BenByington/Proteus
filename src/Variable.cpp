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

#include "Variable.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Solenoid.h"
#include "Field2.h"
#include "Logs/Log.h"

using namespace std;

Variable::MathNode::MathNode() : opName("Empty operator"), opSym("Empty Symbol")
{
}

string Variable::MathNode::executeText()
{
    this->updateLabel();
    
    string ret = getID() + " = " + getDependString() + " -- [ " + getLabel() + " ]";
    
    return ret;
}

Variable::UnaryOp::UnaryOp(GNode* a)
{
    r = a;
}

Variable::BinaryOp::BinaryOp(GNode* a, GNode* b)
{
    l = a;
    r = b;
}

string Variable::UnaryOp::getDependString()
{
    return opName + ": " + r->getID();
}

string Variable::BinaryOp::getDependString()
{
    return opName + ": " + l->getID() + " " + r->getID();
}

void Variable::UnaryOp::updateLabel()
{
    if(presedence < r->getPresedence())
    {
        setLabel(opSym + " (" + r->getLabel() + ")");
    }
    else
    {
        setLabel(opSym + " " + r->getLabel());
    }
}

void Variable::BinaryOp::updateLabel()
{
    string label = "";
    if(presedence < l->getPresedence())
        label += "(" + l->getLabel() + ") ";
    else
        label += l->getLabel() + " ";
    label += opSym;
    if(presedence < r->getPresedence())
        label += " (" + r->getLabel() + ")";
    else
        label += " " + r->getLabel();
    
    setLabel(label);
}

template<class T>
Variable::OperatorTier1<T>::OperatorTier1(shared_ptr<T> left, operation op) : Variable::UnaryOp(left->op)
{
    static_assert(std::is_base_of<Variable,T>::value,"Must use Variable class\n");
    
    presedence = 1;
    switch(op)
    {
        case decompose:
            opName = "Decompose Solenoidal";
            opSym = " : ";
            break;
        case decomposeCurl:
            opName = "Decompose Curled Solenoidal";
            opSym = " : ";
            break;
        case recompose:
            opName = "Recompose Solenoidal";
            opSym = " : ";
    }
}

template <class T>
Variable::OperatorTier2<T>::OperatorTier2(shared_ptr<T> left, operation op) : Variable::UnaryOp(left->op)
{
    static_assert(std::is_base_of<Variable,T>::value,"Must use Variable class\n");
    
    presedence = 2;
    switch(op)
    {
        case div:
            opName = "Divergence";
            opSym = "\u2207 \u00B7";
            break;
        case grad:
            opName = "Gradient";
            opSym = "\u2207";
            break;
        case curl:
            opName = "Curl";
            opSym = "\u2207 \u2A2F";
            break;
        case laplace:
            opName = "Laplacian";
            opSym = "\u2207^2";
    }
}

template<class T1, class T2>
Variable::OperatorTier3<T1,T2>::OperatorTier3(shared_ptr<T1> left, shared_ptr<T2> right, operation op) 
         : Variable::BinaryOp(left->op, right->op)
{
    
    static_assert(std::is_base_of<Variable,T1>::value,"Must use Variable class\n");
    static_assert(std::is_base_of<Variable,T2>::value,"Must use Variable class\n");
    presedence = 3;
    switch(op)
    {
        case mul:
            opName = "Multiply";
            opSym = "*";
            break;
        case divide:
            opName = "Divide";
            opSym = "/";
            break;
        case dot:
            opName = "Dot";
            opSym = "\u00B7";
            
    }
}

template<class T1, class T2>
Variable::OperatorTier4<T1,T2>::OperatorTier4(shared_ptr<T1> left, shared_ptr<T2> right, operation op)
         : Variable::BinaryOp(left->op, right->op)
{
    static_assert(std::is_base_of<Variable,T1>::value,"Must use Variable class\n");
    static_assert(std::is_base_of<Variable,T2>::value,"Must use Variable class\n");
    
    presedence = 4;
    switch(op)
    {
        case cross:
            opName = "Cross";
            opSym = "\u2A2F";
            break;
        case outer:
            opName = "Outer Product";
            opSym = "\u2297";
    }
}

template<class T1, class T2>
Variable::OperatorTier5<T1,T2>::OperatorTier5(shared_ptr<T1> left,  shared_ptr<T2> right, operation op)
         : Variable::BinaryOp(left->op, right->op)
{
    static_assert(std::is_base_of<Variable,T1>::value,"Must use Variable class\n");
    static_assert(std::is_base_of<Variable,T2>::value,"Must use Variable class\n");
    
    presedence = 5;
    switch(op)
    {
        case add:
            opName = "Add";
            opSym = "+";
            break;
        case sub:
            opName = "Subtract";
            opSym = "-";
    }
}

template <class T>
void Variable::OperatorTier1<T>::execute()
{
    error("execute method not implemented!");
}

template <class T>
void Variable::OperatorTier2<T>::execute()
{
    error("execute method not implemented!");
}

template <class T1, class T2>
void Variable::OperatorTier3<T1,T2>::execute()
{
    error("execute method not implemented!");
}

template <class T1, class T2>
void Variable::OperatorTier4<T1,T2>::execute()
{
    error("execute method not implemented!");
}

template <class T1, class T2>
void Variable::OperatorTier5<T1,T2>::execute()
{
    error("execute method not implemented!");
}

template class Variable::OperatorTier1<Vector>;
template class Variable::OperatorTier1<Tensor>;
template class Variable::OperatorTier1<Solenoid>;
template class Variable::OperatorTier2<Field>;
template class Variable::OperatorTier2<Vector>;
template class Variable::OperatorTier2<Tensor>;
template class Variable::OperatorTier3<Scalar,Vector>;
template class Variable::OperatorTier3<Scalar,Field>;
template class Variable::OperatorTier3<Scalar,Tensor>;
template class Variable::OperatorTier3<Field,Scalar>;
template class Variable::OperatorTier3<Field,Field>;
template class Variable::OperatorTier3<Field,Vector>;
template class Variable::OperatorTier3<Field,Tensor>;
template class Variable::OperatorTier3<Vector,Vector>;
template class Variable::OperatorTier4<Vector,Vector>;
template class Variable::OperatorTier5<Field,Field>;
template class Variable::OperatorTier5<Vector,Vector>;
template class Variable::OperatorTier5<Tensor,Tensor>;