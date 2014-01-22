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
    ScalarMul * node = new ScalarMul(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    ScalarDiv * node = new ScalarDiv(fact, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::add(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldAdd * node = new FieldAdd(r, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::subtract(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldSub * node = new FieldSub(r, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::multiply(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldMul * node = new FieldMul(r, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    FieldDiv * node = new FieldDiv(r, shared_from_this());
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field::ScalarMul::ScalarMul(std::shared_ptr<Scalar> s, std::shared_ptr<Field> f)
        : OperatorTier3(f->op, s->op, mul) {}

Field::ScalarDiv::ScalarDiv(std::shared_ptr<Scalar> s, std::shared_ptr<Field> f)
        : OperatorTier3(f->op, s->op, divide) {}

Field::FieldMul::FieldMul(std::shared_ptr<Field> f1, std::shared_ptr<Field> f2)
        : OperatorTier3(f1->op, f2->op, mul) {}

Field::FieldDiv::FieldDiv(std::shared_ptr<Field> f1, std::shared_ptr<Field> f2)
        : OperatorTier3(f1->op, f2->op, divide) {}

Field::FieldAdd::FieldAdd(std::shared_ptr<Field> f1, std::shared_ptr<Field> f2)
        : OperatorTier5(f1->op, f2->op, add) {}

Field::FieldSub::FieldSub(std::shared_ptr<Field> f1, std::shared_ptr<Field> f2)
        : OperatorTier5(f1->op, f2->op, sub) {}