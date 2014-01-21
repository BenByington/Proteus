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

#include "ExecutionGraph.h"
#include "Vector.h"

#include "Logs/Log.h"

#include <iostream>
#include <list>
using namespace std; 

ExecutionGraph::ExecutionGraph()
{
    linearized = false;
}

void ExecutionGraph::registerHead(GNode* h)
{
    this->ready.push(h);
}

void ExecutionGraph::linearize()
{
    if(!order.empty() || linearized == true)
    {
        error("Linearizing the same graph twice!!!");
    }
    
    list<GNode*> free;
    int count = 0;
    while(!ready.empty())
    {
        GNode * curr = ready.top();
        ready.pop();
        
        free = curr->resolveDependency();
        while(!free.empty())
        {
            ready.push(free.back());
            free.pop_back();
        }
        
        curr->setNum(count);
        order.push_back(curr);
        count++;
        //cerr << "Testing: " << curr->executeText() << endl;;
    }
    
    linearized = true;
    
}

void ExecutionGraph::execute()
{
    if(!linearized)
        linearize();
    
    for(vector<GNode*>::iterator i = order.begin(); i != order.end(); i++)
    {
        GNode * curr = (*i);
        cerr << curr->executeText() << endl;
    }
}
