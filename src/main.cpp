#include <iostream>
#include <chrono>

#include "eigen3/Eigen/Dense"
#include "model.hpp"
#include "diffop.hpp"
#include "equation.hpp"

using namespace std;

void myPrint(Eigen::MatrixXd vec, int Nx, int Ny)
{
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            cout << vec(i * Ny + j) << " ";
        }
        cout << endl;
    }
}

int main()
{
    ModelPML model;
    model.ReadJson();

    SourceData source;
    source.ReadJson();

    Acoustic eq(model);
    // eq.Solve(model, source, 8);
    eq.ForwardModelling(model, source);
    eq.WriteData();

    // test openmp
    // model.Nx_pml = 10001;
    // model.Ny_pml = 20001;

    // Eigen::VectorXd vec1(model.Nx_pml * model.Ny_pml);
    // Eigen::VectorXd vec2(model.Nx_pml * model.Ny_pml);
    // Eigen::VectorXd res(model.Nx_pml * model.Ny_pml);
    // DiffOp diff(model);
    // auto time1 = std::chrono::high_resolution_clock::now();
    // res = diff.DyBackward(vec1);
    // auto time2 = std::chrono::high_resolution_clock::now();
    // res = diff.DyBackward_p(vec2);
    // auto time3 = std::chrono::high_resolution_clock::now();

    // auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1);
    // auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2);
    // std::cout << "代码块执行时间: " << duration1.count() << " 毫秒" << std::endl;
    // std::cout << "代码块执行时间: " << duration2.count() << " 毫秒" << std::endl;
}