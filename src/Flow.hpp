#ifndef SRC_FLOW_HPP_
#define SRC_FLOW_HPP_

#include <memory>
#include <string>

#include "LBModel.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"

class Flow final {
 private:
    std::unique_ptr<LBModel> _lbmodel;
    std::unique_ptr<Domain> _domain;
    std::unique_ptr<Boundary> _boundary;
    uint32_t _num_timesteps;

    void init_flow_couette_2d();
    void init_flow_couette_3d();
    void init_flow_poiseuille_2d();
    void init_flow_poiseuille_3d();

 public:
    explicit Flow(std::string flow_type);
    Flow(Flow&) = delete;
    Flow& operator=(Flow&) = delete;
    ~Flow();
    const Domain *get_domain() const;
    const Boundary *get_boundary() const;
    const LBModel *get_lbmodel() const;
    size_t get_num_timesteps() const;
};

#endif  // SRC_FLOW_HPP_
