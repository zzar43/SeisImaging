#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>

#include "json.hpp"

using json = nlohmann::json;

class BasicInfo
{
public:
    int Nx;
    int Ny;
    int Nt;
    double dx;
    double dy;
    double dt;

    BasicInfo() : Nx(0), Ny(0), Nt(0), dx(0.), dy(0.), dt(0.) {}
    ~BasicInfo() {}
};

class Source : virtual public BasicInfo
{
public:
    int source_num;
    Eigen::MatrixXi source_position;
    Eigen::MatrixXd source_fn;
};

class Recevier : virtual public BasicInfo
{
public:
    int receiver_num;
    Eigen::MatrixXi receiver_position;
};

class Model : public Source, public Recevier
{
public:
    Eigen::MatrixXd c;
    Eigen::MatrixXd rho;
};

class ModelPML : public Model
{
public:
    int Nx_pml;
    int Ny_pml;
    Eigen::MatrixXd c_pml;
    Eigen::MatrixXd rho_pml;
    Eigen::MatrixXd sigma_x;
    Eigen::MatrixXd sigma_y;
    Eigen::MatrixXi source_position_pml;
    Eigen::MatrixXi receiver_position_pml;

    ModelPML() : Nx_pml(0), Ny_pml(0) {}
    ~ModelPML() {}

    void ReadJson();

private:
    void ReadArray2D(Eigen::MatrixXd &var, json data, std::string var_name, int dim1, int dim2);
    void ReadArray2D(Eigen::MatrixXi &var, json data, std::string var_name, int dim1, int dim2);
};