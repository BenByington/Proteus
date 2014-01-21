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

#ifndef EXECUTIONGRAPH_H
#define	EXECUTIONGRAPH_H

#include "GNode.h"

#include <vector>
#include <queue>

class ExecutionGraph final
{
public:
    ExecutionGraph();
    ~ExecutionGraph(){};
    
    void execute();
    void registerHead(GNode * h);
private:
    void linearize();
    bool linearized;
    
    std::vector<GNode*> order;
    std::priority_queue<GNode*> ready;
};

#endif	/* EXECUTIONGRAPH_H */

