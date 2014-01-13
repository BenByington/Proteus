#include "Variable.h"
#include "Scalar.h"

Variable::Variable()
{
    op = 0;
}

Variable * Variable::operator *(Scalar * fact)
{
    Variable * ret = createVariable();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Variable * Variable::operator /(Scalar * fact)
{
    Variable * ret = createVariable();
    ScalarFactor * node = new ScalarFactor(fact, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(fact->op);
    
    return ret;
}

Variable * Variable::operator +(Variable * r)
{
    Variable * ret = createVariable();
    VariableArithmetic * node = new VariableArithmetic(r, this);
    node->op = node->add;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Variable * Variable::operator -(Variable * r)
{
    Variable * ret = createVariable();
    VariableArithmetic * node = new VariableArithmetic(r, this);
    node->op = node->sub;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Variable * Variable::operator *(Variable * r)
{
    Variable * ret = createVariable();
    VariableArithmetic * node = new VariableArithmetic(r, this);
    node->op = node->mul;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

Variable * Variable::operator /(Variable * r)
{
    Variable * ret = createVariable();
    VariableArithmetic * node = new VariableArithmetic(r, this);
    node->op = node->divide;
    ret->op = node;
    
    node->addDependency(this->op);
    node->addDependency(r->op);
    
    return ret;
}

void Variable::ScalarFactor::execute()
{
    
}

void Variable::VariableArithmetic::execute()
{
    
}