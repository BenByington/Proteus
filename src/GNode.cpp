#include "GNode.h"

#include "Logs/Log.h"

#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

GNode::GNode() : label("")
{
    presedence = 0;
}

void GNode::setLabel(std::string s)
{
    label = s;
}

void GNode::addDependency(GNode * p)
{
    parents.push_back(p);
    p->children.push_back(this);
}

string GNode::getLabel()
{
    return label;
}

string GNode::getID()
{
    stringstream ss;
    ss << "N" << this->myNum;
    return ss.str();
}

void GNode::setNum(int n)
{
    myNum = n;
}

list<GNode*> GNode::resolveDependency()
{
    if(!parents.empty())
    {
        error("Resolving GNode when parent nodes have not yet been resolved!");
    }
    
    list<GNode*> ret;
    //get rid of duplicates
    children.sort();
    children.unique();
    
    for(list<GNode*>::iterator i = children.begin(); i != children.end(); i++)
    {
        GNode * child = *i;
        
        int c1 = child->parents.size();
        child->parents.remove(this);
        int c2 = child->parents.size();
        
        if(c1 != c2-1)
        {
            error("We did not find ourself as a dependency in one of our children!")
        }
        
        if(c2 == 0)
            ret.push_back(child);
    }
    
    return ret;
}

string GNode::executeText()
{
    //this->updateLabel();
   
    string ret = getID() + " " + getLabel(); //+ " = " + getDependString() + " -- [ " + getLabel() + " ]";
    
    return ret;
}