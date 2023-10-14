#include "equation.hpp"

MatrixXd Equation2D::dx_forward(MatrixXd U, double h)
{
    MatrixXd Ux(U.rows(), U.cols());
    for (int j = 1; j < U.cols() - 1; j++)
    {
        for (int i = 1; i < U.rows() - 1; i++)
        {
            Ux(i, j) = (U(i + 1, j) - U(i, j)) / h;
        }
    }
    return Ux;
}

MatrixXd Equation2D::dx_backward(MatrixXd U, double h)
{
    MatrixXd Ux(U.rows(), U.cols());
    for (int j = 1; j < U.cols() - 1; j++)
    {
        for (int i = 1; i < U.rows() - 1; i++)
        {
            Ux(i, j) = (U(i, j) - U(i - 1, j)) / h;
        }
    }
    return Ux;
}

MatrixXd Equation2D::dy_forward(MatrixXd U, double h)
{
    MatrixXd Ux(U.rows(), U.cols());
    for (int j = 1; j < U.cols() - 1; j++)
    {
        for (int i = 1; i < U.rows() - 1; i++)
        {
            Ux(i, j) = (U(i, j + 1) - U(i, j)) / h;
        }
    }
    return Ux;
}

MatrixXd Equation2D::dy_backward(MatrixXd U, double h)
{
    MatrixXd Ux(U.rows(), U.cols());
    for (int j = 1; j < U.cols() - 1; j++)
    {
        for (int i = 1; i < U.rows() - 1; i++)
        {
            Ux(i, j) = (U(i, j) - U(i, j - 1)) / h;
        }
    }
    return Ux;
}


void AcousticWaveEq2D::Solve(int source_idx)
{
    // initialization
    a = MatrixXd::Zero(Nx_pml, Ny_pml);
    b = MatrixXd::Zero(Nx_pml, Ny_pml);
    u1 = MatrixXd::Zero(Nx_pml, Ny_pml);
    u2 = MatrixXd::Zero(Nx_pml, Ny_pml);
    vx1 = MatrixXd::Zero(Nx_pml, Ny_pml);
    vx2 = MatrixXd::Zero(Nx_pml, Ny_pml);
    vy1 = MatrixXd::Zero(Nx_pml, Ny_pml);
    vy2 = MatrixXd::Zero(Nx_pml, Ny_pml);
    phi1 = MatrixXd::Zero(Nx_pml, Ny_pml);
    phi2 = MatrixXd::Zero(Nx_pml, Ny_pml);
    psi1 = MatrixXd::Zero(Nx_pml, Ny_pml);
    psi2 = MatrixXd::Zero(Nx_pml, Ny_pml);
    dxvx1_f = MatrixXd::Zero(Nx_pml, Ny_pml);
    dyvy1_f = MatrixXd::Zero(Nx_pml, Ny_pml);
    dxvx1_b = MatrixXd::Zero(Nx_pml, Ny_pml);
    dyvy1_b = MatrixXd::Zero(Nx_pml, Ny_pml);
    dxu2 = MatrixXd::Zero(Nx_pml, Ny_pml);
    dyu2 = MatrixXd::Zero(Nx_pml, Ny_pml);

    MatrixXd part1 = MatrixXd::Zero(Nx_pml, Ny_pml);
    MatrixXd part2 = MatrixXd::Zero(Nx_pml, Ny_pml);

    a = rho_pml.cwiseInverse();
    b = rho_pml.cwiseProduct(c_pml.cwiseProduct(c_pml));

    for (int iter = 0; iter < Nt; iter++)
    {
        dxvx1_f = dx_forward(vx1, dx);
        dyvy1_f = dy_forward(vy1, dy);
        dxvx1_b = dx_backward(vx1, dx);
        dyvy1_b = dy_backward(vy1, dy);

        part1 = -1 * u1.cwiseProduct(sigma_x + sigma_y);
        part2 = b.cwiseProduct(dxvx1_f + dyvy1_f);
        u2 = u1 + dt * (part1 + part2 + phi1 + psi1);

        // update source
        u2(source_position_pml(source_idx, 0), source_position_pml(source_idx, 1)) += b(source_position_pml(source_idx, 0), source_position_pml(source_idx, 1)) * source_fn(source_idx, iter) * dt;

        dxu2 = dx_backward(u2, dx);
        dyu2 = dy_backward(u2, dy);

        vx2 = vx1 + dt * (a.cwiseProduct(dxu2) - sigma_x.cwiseProduct(vx1));
        vy2 = vy1 + dt * (a.cwiseProduct(dyu2) - sigma_y.cwiseProduct(vy1));
        phi2 = phi1 + dt * b.cwiseProduct(sigma_y.cwiseProduct(dxvx1_f));
        psi2 = psi1 + dt * b.cwiseProduct(sigma_x.cwiseProduct(dyvy1_f));

        // update
        u1 = u2;
        vx1 = vx2;
        vy1 = vy2;
        phi1 = phi2;
        psi1 = psi2;
    }
}

void AcousticWaveEq2D::WriteJson()
{
    json j_u2 = AcousticWaveEq2D::EigenToJson(u2);

    json output = {
        {"u2", j_u2}
    };

    std::ofstream o("./data/temp_cpp.json");
    o << std::setw(4) << output << std::endl;
}

json AcousticWaveEq2D::EigenToJson(MatrixXd A)
{
    std::vector<double> v(A.data(), A.data() + A.rows() * A.cols());
    json j_vec(v);
    return j_vec;
}