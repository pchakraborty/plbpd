#ifndef FLOW_HPP
#define FLOW_HPP

#include <memory>
#include "Domain.hpp"
#include "Boundary.hpp"

class Flow final{

private:

    std::unique_ptr<Domain> _domain;
    std::unique_ptr<Boundary> _boundary;

    void initCouetteFlow();
   
public:

    Flow(std::string flowType);
    ~Flow();
    const Domain *getFlowDomain() const;
    const Boundary *getFlowBoundary() const;
    
};

#endif
