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

#include <vector>
#include <string>

class GNode
{
protected:
    
std::vector<GNode*> parents; 
std::vector<GNode*> children;

public:
    virtual ~GNode(){}
	
    virtual void execute() = 0;
    virtual std::string executeText() = 0;
    
    void addDependency(GNode * p);
    std::string getName();
    
    int myNum;
protected:
    static int numNodes;
};

#endif