#pragma once

#include <string>
#include <eigen3/Eigen/Dense>

#include "model.hpp"

using namespace Eigen;

class Equation2D : public ModelPML
{
public:
    MatrixXd dx_forward(MatrixXd U, double h);
    MatrixXd dx_backward(MatrixXd U, double h);
    MatrixXd dy_forward(MatrixXd U, double h);
    MatrixXd dy_backward(MatrixXd U, double h);
};

class AcousticWaveEq2D : public Equation2D
{
public:
    MatrixXd U, data;

    AcousticWaveEq2D() {}
    ~AcousticWaveEq2D() {}
    
    void Solve(int source_idx=0);

    void WriteJson();
    json EigenToJson(MatrixXd A);
private:
    MatrixXd a, b, u1, u2, vx1, vx2, vy1, vy2, phi1, phi2, psi1, psi2, dxvx1_f, dyvy1_f, dxvx1_b, dyvy1_b, dxu2, dyu2;

};