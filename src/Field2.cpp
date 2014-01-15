#include "Field2.h"
#include "Scalar.h"

Field::Field()
{
    op = 0;
}

Field * Field::operator *(Scalar * fact)
{
    Field * ret = createField();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Field * Field::operator /(Scalar * fact)
{
    Field * ret = createField();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Field * Field::operator +(Field * r)
{
    Field * ret = createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Field::operator -(Field * r)
{
    Field * ret = createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Field::operator *(Field * r)
{
    Field * ret = createField();
    FieldArithmetic * node = new FieldArithmetic(r, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Field::operator /(Field * r)
{
    Field * ret = createField();
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

Field::FieldArithmetic::FieldArithmetic(Field* p1, Field* p2)
{
    this->p1 = p1;
    this->p2 = p2;
}

void Field::FieldArithmetic::execute()
{
    
}