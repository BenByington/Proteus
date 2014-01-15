#include "Vector.h"
#include "Scalar.h"

Vector::Vector()
{
    op = 0;
}

Vector * Vector::operator +(Vector * r)
{
    Vector * ret = createVector();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Vector * Vector::operator -(Vector * r)
{
    Vector * ret = createVector();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Vector * Vector::cross(Vector * r)
{
    Vector * ret = createVector();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->cross;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Field * Vector::dot(Vector* r)
{
    Field * ret = createField();
    VectorArithmetic * node = new VectorArithmetic(this, r);
    node->op = node->dot;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
    
}

Vector::VectorArithmetic::VectorArithmetic(Vector* v1, Vector* v2)
{
    this->p1 = v1;
    this->p2 = v2;
}

void Vector::VectorArithmetic::execute()
{
    
}
