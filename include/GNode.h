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

#ifndef	GNODE_H
#define GNODE_H

#include <list>
#include <string>

class GNode
{
protected:
    
std::list<GNode*> parents; 
std::list<GNode*> children;

public:
    virtual ~GNode(){}
	
    virtual void execute() = 0;
    std::string executeText();
    
    void addDependency(GNode * p);
    std::list<GNode*> resolveDependency();
    std::string getLabel();
    std::string getID();
    virtual std::string getDependString() = 0;
    
    void setNum(int n);
    void setLabel(std::string s);
    
protected:
    GNode();
    int myNum;
    std::string label;
};

#endif