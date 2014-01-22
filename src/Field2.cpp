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
    OperatorTier3<Scalar,Field> * node = new OperatorTier3<Scalar,Field>(fact,shared_from_this(),OperatorTier3<Scalar,Field>::mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Scalar> fact)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier3<Field,Scalar> * node = new OperatorTier3<Field,Scalar>(shared_from_this(),fact,OperatorTier3<Field,Scalar>::divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

shared_ptr<Field> Field::add(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier5<Field,Field> * node = new OperatorTier5<Field,Field>(shared_from_this(),r,OperatorTier5<Field,Field>::add);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::subtract(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier5<Field,Field> * node = new OperatorTier5<Field,Field>(shared_from_this(),r,OperatorTier5<Field,Field>::sub);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::multiply(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier3<Field,Field> * node = new OperatorTier3<Field,Field>(shared_from_this(),r,OperatorTier3<Field,Field>::mul);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

shared_ptr<Field> Field::divide(shared_ptr<Field> r)
{
    shared_ptr<Field> ret = VariableFactory::createField();
    OperatorTier3<Field,Field> * node = new OperatorTier3<Field,Field>(shared_from_this(),r,OperatorTier3<Field,Field>::divide);
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}
