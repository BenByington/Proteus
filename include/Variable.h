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

class Vector;
class Scalar;

class Variable
{
protected:
    Variable(){op = new nop();}
    
    class nop : public GNode
    {
        friend class Variable;
    public:
        nop(){};
        virtual void execute(){};
        virtual std::string getDependString(){return "";}
    };
    
public:
    virtual ~Variable(){}
    
    GNode * op;

};

#endif