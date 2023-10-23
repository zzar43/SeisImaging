#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Dense"
#include "json.hpp"

using json = nlohmann::json;

// container
class SeisData
{
public:
    int Nx_pml, Ny_pml, Nt;
    Eigen::VectorXd last_wavefield, wavefield, seis_signal, seis_signal_multi;
};

class SourceData
{
public:
    int source_num, Nt;
    // data should be a 2d matrix with source_num by Nt
    Eigen::VectorXd data;
    void ReadJson();
    void ReadVector(json data, Eigen::VectorXd &vec, std::string var_name, int N);
};


class Mesh
{
public:
    int Nx, Ny;
    double dx, dy;
    Eigen::VectorXd x, y;
};

class Time
{
public:
    int Nt;
    double dt;
    Eigen::VectorXd t;
};

class Acquisition
{
public:
    int source_num, receiver_num;
    Eigen::VectorXi source_position, receiver_position;
};

class Model
{
public:
    int Nx, Ny, Nt;
    double dx, dy, dt;
    Eigen::VectorXd x, y, t;

    int source_num, receiver_num;
    Eigen::VectorXi source_position, receiver_position;

    // fields
    Eigen::VectorXd c, rho;

    Model() {}
    Model(const Mesh &mesh, const Time &time, const Acquisition &acquisition);

    void ReadVelocity(){};
    void ReadDensity(){};
};

class ModelPML : public Model
{
public:
    // read data
    json data;

    // PML
    int pml_len;
    double pml_alpha;

    int Nx_pml, Ny_pml;
    Eigen::VectorXi source_position_pml, receiver_position_pml;

    // fields
    Eigen::VectorXd c_pml, rho_pml, sigma_x, sigma_y;

    // SeisData
    SeisData seis_data;

    ModelPML();
    ModelPML(const Mesh &mesh, const Time &time, const Acquisition &acquisition, int pml_len, double pml_alpha);

    void ReadJson();
    void ReadVector(json data, Eigen::VectorXd &vec, std::string var_name, int N);
    void ReadVector(json data, Eigen::VectorXi &vec, std::string var_name, int N);
};