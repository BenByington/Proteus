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
    node->setOp(node->mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    ScalarFactor * node = new ScalarFactor(fact, shared_from_this());
    node->setOp(node->divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::add(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->setOp(node->add);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::subtract(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->setOp(node->sub);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::multiply(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->setOp(node->mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldArithmetic * node = new FieldArithmetic(r, shared_from_this());
    node->setOp(node->divide);
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

string Field::ScalarFactor::getDependString()
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
    
    string ret = opName + string(": ") + this->sParent->op->getID();
    ret += string(" ") + this->vParent->op->getID();
    
    
    return ret;
}

void Field::ScalarFactor::setOp(operations o)
{
    op = o;
    string opName;
    switch(op)
    {
    case mul:
        opName = " * ";
        break;
    case divide:
        opName = " / ";
        break;
    }
    
    label = "(" + vParent->op->getLabel() + opName + sParent->op->getLabel() + ")";
    
}

Field::FieldArithmetic::FieldArithmetic(shared_ptr<Field> p1, shared_ptr<Field> p2)
{
    this->p1 = p1;
    this->p2 = p2;
}

void Field::FieldArithmetic::execute()
{
    
}

string Field::FieldArithmetic::getDependString()
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
    
    string ret = opName + ": " + this->p1->op->getID();
    ret += " " + this->p2->op->getID();
    
    return ret;
}

void Field::FieldArithmetic::setOp(operations o)
{
    op = o;
    string opName;
    switch(op)
    {
    case mul:
        opName = " * ";
        break;
    case divide:
        opName = " / ";
        break;
    case sub:
        opName = " - ";
        break;
    case add:
        opName = " + ";    
    }
    
    label = "(" + p1->op->getLabel() + opName + p2->op->getLabel() + ")";
    
}