#include <iostream>
#include <eigen3/Eigen/Dense>
#include "model.hpp"

using namespace Eigen;
using namespace std;

void AddOne(MatrixXd &A)
{
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 20; j++)
        {
            A(i, j) += 1;
        }
    }
}

int main()
{

    ModelPML model;

    model.ReadJson();

    cout << model.Nx << endl;
    cout << model.Ny << endl;
    cout << model.c_pml << endl;
    cout << model.source_position_pml << endl;
    cout << model.receiver_position_pml << endl;
    cout << model.source_fn << endl;


    // MatrixXd A(10,20);
    // for (int i = 0; i < 10; i++)
    // {
    //     for (int j = 0; j < 20; j++)
    //     {
    //         A(i, j) = i + j;
    //     }
    // }

    // cout << A << endl;

    // AddOne(A);

    // cout << A << endl;
}