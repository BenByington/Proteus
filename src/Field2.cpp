#include "Field2.h"
#include "Scalar.h"
#include "VariableFactory.h"

using namespace std;

Field::Field()
{
}

shared_ptr<Field> Field::multiply(shared_ptr<Scalar> fact)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::add(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::subtract(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::multiply(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field::ScalarFactor::ScalarFactor(shared_ptr<Scalar> s, shared_ptr<Field> v)
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

Field::FieldArithmetic::FieldArithmetic(shared_ptr<Field> p1, shared_ptr<Field> p2)
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