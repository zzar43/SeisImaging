#pragma once

#include "model.hpp"
#include "omp.h"
#include "eigen3/Eigen/Dense"

// differential operator in 2d
class DiffOp
{
private:
    int Nx, Ny;
    double dx, dy;

public:
    DiffOp() {}
    DiffOp(const ModelPML &model) : Nx(model.Nx_pml), Ny(model.Ny_pml), dx(model.dx), dy(model.dy){};

    Eigen::VectorXd DxForward(const Eigen::VectorXd &vec);
    Eigen::VectorXd DxBackward(const Eigen::VectorXd &vec);
    Eigen::VectorXd DyForward(const Eigen::VectorXd &vec);
    Eigen::VectorXd DyBackward(const Eigen::VectorXd &vec);
    Eigen::VectorXd DyBackward_p(const Eigen::VectorXd &vec);
};