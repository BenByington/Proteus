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

Variable::OperatorTier1::OperatorTier1(GNode * left, operation op) : Variable::UnaryOp(left)
{
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

Variable::OperatorTier2::OperatorTier2(GNode * left, operation op) : Variable::UnaryOp(left)
{
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

Variable::OperatorTier3::OperatorTier3(GNode * left, GNode * right, operation op) 
         : Variable::BinaryOp(left, right)
{
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

Variable::OperatorTier4::OperatorTier4(GNode * left, GNode * right, operation op)
         : Variable::BinaryOp(left, right)
{
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

Variable::OperatorTier5::OperatorTier5(GNode * left,  GNode * right, operation op)
         : Variable::BinaryOp(left, right)
{
    
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