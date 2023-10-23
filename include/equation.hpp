#pragma once

// #include <vector>

#include "model.hpp"
#include "diffop.hpp"

#include "eigen3/Eigen/Dense"

class Equation
{
public:
    virtual void Solve(){};
    virtual void ForwardModelling();
    virtual void BackwardOperator();
    virtual void AdjointMethod();
};

class Acoustic
{
private:
    DiffOp diff;

public:
    // this seis_data is for compute adjoint method
    SeisData seis_data;

    Acoustic(){};
    Acoustic(const ModelPML &model) : diff(DiffOp(model)){};
    ~Acoustic(){};

    void Solve(const ModelPML &model, const SourceData &source_data, int source_idx = 0);
    void TimeUpdate(int time_iter, const ModelPML &model, const SourceData &source_data, int source_idx);
    void ForwardModelling(const ModelPML &model, const SourceData &source_data);

    void RecordSeisSignal(int time_iter, const ModelPML &model);

    void WriteData();
    json EigenToJson(const Eigen::VectorXd &vec);

private:
    Eigen::VectorXd a, b, u1, u2, vx1, vx2, vy1, vy2, phi1, phi2, psi1, psi2, dxvx1_f, dyvy1_f, dxvx1_b, dyvy1_b, dxu2, dyu2, part1, part2;
};