#include "GNode.h"

using namespace std;

void GNode::addDependency(GNode * p)
{
    parents.push_back(p);
    p->children->push_back(this);
}