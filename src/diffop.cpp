#include "diffop.hpp"
#include <iostream>

// 2d index
// #define idx(i, j, Ny) (i * Ny + j)

// inline int idxi(int i, int j, int Ny)
// {
//     return i * Ny + j;
// }

Eigen::VectorXd DiffOp::DxForward(const Eigen::VectorXd &vec)
{
    Eigen::VectorXd vec_x = Eigen::VectorXd::Zero(vec.rows());
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < Nx - 1; i++)
    {
        for (int j = 1; j < Ny - 1; j++)
        {
            vec_x(i * Ny + j) = (vec((i + 1) * Ny + j) - vec(i * Ny + j)) / dx;
        }
    }
    return vec_x;
}

Eigen::VectorXd DiffOp::DxBackward(const Eigen::VectorXd &vec)
{
    Eigen::VectorXd vec_x = Eigen::VectorXd::Zero(vec.rows());
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < Nx - 1; i++)
    {
        for (int j = 1; j < Ny - 1; j++)
        {
            vec_x(i * Ny + j) = (vec(i * Ny + j) - vec((i - 1) * Ny + j)) / dx;
        }
    }
    return vec_x;
}

Eigen::VectorXd DiffOp::DyForward(const Eigen::VectorXd &vec)
{
    Eigen::VectorXd vec_y = Eigen::VectorXd::Zero(vec.rows());
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < Nx - 1; i++)
    {
        for (int j = 1; j < Ny - 1; j++)
        {
            vec_y(i * Ny + j) = (vec(i * Ny + j + 1) - vec(i * Ny + j)) / dy;
        }
    }
    return vec_y;
}

Eigen::VectorXd DiffOp::DyBackward(const Eigen::VectorXd &vec)
{
    Eigen::VectorXd vec_y = Eigen::VectorXd::Zero(vec.rows());
    // #pragma omp parallel for collapse(2)
    for (int i = 1; i < Nx - 1; i++)
    {
        for (int j = 1; j < Ny - 1; j++)
        {
            vec_y(i * Ny + j) = (vec(i * Ny + j) - vec(i * Ny + j - 1)) / dy;
        }
    }
    return vec_y;
}

Eigen::VectorXd DiffOp::DyBackward_p(const Eigen::VectorXd &vec)
{
    Eigen::VectorXd vec_y = Eigen::VectorXd::Zero(vec.rows());
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < Nx - 1; i++)
    {
        for (int j = 1; j < Ny - 1; j++)
        {
            vec_y(i * Ny + j) = (vec(i * Ny + j) - vec(i * Ny + j - 1)) / dy;
        }
    }
    return vec_y;
}
