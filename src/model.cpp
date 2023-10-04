#include "model.hpp"

void ModelPML::ReadJson()
{
    std::ifstream f("sample.json");
    json data = json::parse(f);

    Nx = data["Nx"];
    Ny = data["Ny"];
    Nt = data["Nt"];
    dx = data["dx"];
    dy = data["dy"];
    dt = data["dt"];
    source_num = data["source_num"];
    receiver_num = data["receiver_num"];
    Nx_pml = data["Nx_pml"];
    Ny_pml = data["Ny_pml"];

    c_pml = Eigen::MatrixXd(Nx_pml, Ny_pml);
    ReadArray2D(c_pml, data, "c_pml", Nx_pml, Ny_pml);
    rho_pml = Eigen::MatrixXd(Nx_pml, Ny_pml);
    ReadArray2D(rho_pml, data, "rho_pml", Nx_pml, Ny_pml);

    source_position_pml = Eigen::MatrixXi(source_num, 2);
    ReadArray2D(source_position_pml, data, "source_position_pml", source_num, 2);
    source_fn = Eigen::MatrixXd(source_num, Nt);
    ReadArray2D(source_fn, data, "source_fn", source_num, Nt);
    receiver_position_pml = Eigen::MatrixXi(receiver_num, 2);
    ReadArray2D(receiver_position_pml, data, "receiver_position_pml", receiver_num, 2);

    sigma_x = Eigen::MatrixXd(Nx_pml, Ny_pml);
    ReadArray2D(sigma_x, data, "sigma_x", Nx_pml, Ny_pml);
    sigma_y = Eigen::MatrixXd(Nx_pml, Ny_pml);
    ReadArray2D(sigma_y, data, "sigma_y", Nx_pml, Ny_pml);
}

void ModelPML::ReadArray2D(Eigen::MatrixXd &var, json data, std::string var_name, int dim1, int dim2)
{
    for (int i = 0; i < dim1; i++)
    {
        for (int j = 0; j < dim2; j++)
        {
            var(i, j) = data[var_name][i][j];
        }
    }
}

void ModelPML::ReadArray2D(Eigen::MatrixXi &var, json data, std::string var_name, int dim1, int dim2)
{
    for (int i = 0; i < dim1; i++)
    {
        for (int j = 0; j < dim2; j++)
        {
            var(i, j) = data[var_name][i][j];
        }
    }
}