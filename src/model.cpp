#include "model.hpp"

#include <iostream>
using namespace std;

// container
void SourceData::ReadJson()
{
    std::ifstream f("./data/source.json");
    json j_data = json::parse(f);

    source_num = j_data["source_num"];
    Nt = j_data["Nt"];
    data = Eigen::VectorXd(source_num * Nt);
    ReadVector(j_data, data, "data", source_num * Nt);
};

void SourceData::ReadVector(json data, Eigen::VectorXd &vec, std::string var_name, int N)
{
    for (int i = 0; i < N; i++)
    {
        vec(i) = data[var_name][i];
    }
}

Model::Model(const Mesh &mesh, const Time &time, const Acquisition &acquisition)
{
    Nx = mesh.Nx;
    Ny = mesh.Ny;
    dx = mesh.dx;
    dy = mesh.dy;

    Nt = time.Nt;
    dt = time.dt;

    source_num = acquisition.source_num;
    receiver_num = acquisition.receiver_num;
    source_position = acquisition.source_position;
    receiver_position = acquisition.receiver_position;
};

ModelPML::ModelPML(const Mesh &mesh, const Time &time, const Acquisition &acquisition, int pml_len, double pml_alpha)
{
    Nx = mesh.Nx;
    Ny = mesh.Ny;
    dx = mesh.dx;
    dy = mesh.dy;

    Nt = time.Nt;
    dt = time.dt;

    source_num = acquisition.source_num;
    receiver_num = acquisition.receiver_num;
    source_position = acquisition.source_position;
    receiver_position = acquisition.receiver_position;
};

ModelPML::ModelPML() {}

void ModelPML::ReadJson()
{
    std::ifstream f("./data/data.json");
    json data = json::parse(f);

    Nx = data["Nx"];
    Ny = data["Ny"];
    Nx_pml = data["Nx_pml"];
    Ny_pml = data["Ny_pml"];
    Nt = data["Nt"];
    dx = data["dx"];
    dy = data["dy"];
    dt = data["dt"];
    source_num = data["source_num"];
    receiver_num = data["receiver_num"];

    x = Eigen::VectorXd(Nx);
    y = Eigen::VectorXd(Ny);
    t = Eigen::VectorXd(Nt);

    c_pml = Eigen::VectorXd(Nx_pml * Ny_pml);
    rho_pml = Eigen::VectorXd(Nx_pml * Ny_pml);
    sigma_x = Eigen::VectorXd(Nx_pml * Ny_pml);
    sigma_y = Eigen::VectorXd(Nx_pml * Ny_pml);
    source_position_pml = Eigen::VectorXi(2 * source_num);
    receiver_position_pml = Eigen::VectorXi(2 * receiver_num);

    ReadVector(data, x, "x", Nx);
    ReadVector(data, y, "y", Ny);
    ReadVector(data, t, "t", Nt);

    ReadVector(data, c_pml, "c_pml", Nx_pml * Ny_pml);
    ReadVector(data, rho_pml, "rho_pml", Nx_pml * Ny_pml);
    ReadVector(data, sigma_x, "sigma_x", Nx_pml * Ny_pml);
    ReadVector(data, sigma_y, "sigma_y", Nx_pml * Ny_pml);
    ReadVector(data, source_position_pml, "source_position_pml", 2 * source_num);
    ReadVector(data, receiver_position_pml, "receiver_position_pml", 2 * receiver_num);
}

void ModelPML::ReadVector(json data, Eigen::VectorXd &vec, std::string var_name, int N)
{
    for (int i = 0; i < N; i++)
    {
        vec(i) = data[var_name][i];
    }
}

void ModelPML::ReadVector(json data, Eigen::VectorXi &vec, std::string var_name, int N)
{
    for (int i = 0; i < N; i++)
    {
        vec(i) = data[var_name][i];
    }
}