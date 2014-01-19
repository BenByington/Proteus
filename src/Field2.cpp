#include "Field2.h"
#include "Scalar.h"
#include "VariableFactory.h"

using namespace std;

Field::Field()
{
}

Field * Field::multiply(Scalar * fact)
{
    Field * ret = VariableFactory::createField();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Field * Field::divide(Scalar * fact)
{
    Field * ret = VariableFactory::createField();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Field * Field::add(Field * r)
{
    Field * ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Field::subtract(Field * r)
{
    Field * ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Field::multiply(Field * r)
{
    Field * ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Field::divide(Field * r)
{
    Field * ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field::ScalarFactor::ScalarFactor(Scalar* s, Field* v)
{
    sParent = s;
    vParent = v;
}

void Field::ScalarFactor::execute()
{
    
}

string Field::ScalarFactor::executeText()
{
    string opName;
    switch(op)
    {
    case mul:
        opName = "Multiply";
        break;
    case divide:
        opName = "Divide";
        break;
    }
    
    string ret = getName() + string(" = "); 
    ret += opName + string(": ") + this->sParent->op->getName();
    ret += string(" ") + this->vParent->op->getName();
    
    
    return ret;
}

Field::FieldArithmetic::FieldArithmetic(Field* p1, Field* p2)
{
    this->p1 = p1;
    this->p2 = p2;
}

void Field::FieldArithmetic::execute()
{
    
}

string Field::FieldArithmetic::executeText()
{
    string opName;
    switch(op)
    {
    case mul:
        opName = "Multiply";
        break;
    case divide:
        opName = "Divide";
        break;
    case sub:
        opName = "Subtract";
        break;
    case add:
        opName = "Add";    
    }
    
    string ret = getName() + " = "; 
    ret += opName + ": " + this->p1->op->getName();
    ret += " " + this->p2->op->getName();
    
    return ret;
}