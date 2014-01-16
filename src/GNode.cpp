#include "GNode.h"

#include <sstream>

using namespace std;

void GNode::addDependency(GNode * p)
{
    parents.push_back(p);
    p->children.push_back(this);
}

string GNode::getName()
{
    stringstream ss;
    ss << "N" << this->myNum;
    return ss.str();
}

int GNode::numNodes = 0;