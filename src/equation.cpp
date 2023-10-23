#include "equation.hpp"

void Acoustic::Solve(const ModelPML &model, const SourceData &source_data, int source_idx)
{
    // init
    a = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    b = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    u1 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    u2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    vx1 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    vx2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    vy1 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    vy2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    phi1 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    phi2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    psi1 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    psi2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    dxvx1_f = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    dyvy1_f = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    dxvx1_b = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    dyvy1_b = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    dxu2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    dyu2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);

    part1 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);
    part2 = Eigen::VectorXd::Zero(model.Nx_pml * model.Ny_pml);

    a = model.rho_pml.cwiseInverse();
    b = model.rho_pml.cwiseProduct(model.c_pml.cwiseProduct(model.c_pml));

    // record
    seis_data.seis_signal = Eigen::VectorXd::Zero(model.receiver_num * model.Nt);

    for (int iter = 0; iter < model.Nt; iter++)
    {
        TimeUpdate(iter, model, source_data, source_idx);
        RecordSeisSignal(iter, model);
    }
    seis_data.last_wavefield = u2;
};

void Acoustic::TimeUpdate(int time_iter, const ModelPML &model, const SourceData &source_data, int source_idx)
{
    int s_idx_1 = model.source_position_pml(source_idx * 2 + 0);
    int s_idx_2 = model.source_position_pml(source_idx * 2 + 1);

    int Ny = model.Ny_pml;

    dxvx1_f = diff.DxForward(vx1);
    dyvy1_f = diff.DyForward(vy1);
    dxvx1_b = diff.DxBackward(vx1);
    dyvy1_b = diff.DyBackward(vy1);

    part1 = -1 * u1.cwiseProduct(model.sigma_x + model.sigma_y);
    part2 = b.cwiseProduct(dxvx1_f + dyvy1_f);
    u2 = u1 + model.dt * (part1 + part2 + phi1 + psi1);

    // update source
    u2(s_idx_1 * Ny + s_idx_2) += b(s_idx_1 * Ny + s_idx_2) * source_data.data(source_idx * model.Nt + time_iter) * model.dt;

    dxu2 = diff.DxBackward(u2);
    dyu2 = diff.DyBackward(u2);

    vx2 = vx1 + model.dt * (a.cwiseProduct(dxu2) - model.sigma_x.cwiseProduct(vx1));
    vy2 = vy1 + model.dt * (a.cwiseProduct(dyu2) - model.sigma_y.cwiseProduct(vy1));
    phi2 = phi1 + model.dt * b.cwiseProduct(model.sigma_y.cwiseProduct(dxvx1_f));
    psi2 = psi1 + model.dt * b.cwiseProduct(model.sigma_x.cwiseProduct(dyvy1_f));

    // update
    u1 = u2;
    vx1 = vx2;
    vy1 = vy2;
    phi1 = phi2;
    psi1 = psi2;
};

void Acoustic::ForwardModelling(const ModelPML &model, const SourceData &source_data)
{
    int N = model.receiver_num * model.Nt;
    seis_data.seis_signal_multi = Eigen::VectorXd::Zero(model.source_num * N);

    for (int source_idx = 0; source_idx < model.source_num; source_idx++)
    {
        Solve(model, source_data, source_idx);
        for (int i = 0; i < N; i++)
        {
            seis_data.seis_signal_multi(source_idx * N + i) = seis_data.seis_signal(i);
        }
    }
}

void Acoustic::RecordSeisSignal(int time_iter, const ModelPML &model)
{
    int Ny = model.Ny_pml;
    for (int receiver_idx = 0; receiver_idx < model.receiver_num; receiver_idx++)
    {
        int r_idx_1 = model.receiver_position_pml(receiver_idx * 2 + 0);
        int r_idx_2 = model.receiver_position_pml(receiver_idx * 2 + 1);
        seis_data.seis_signal(receiver_idx * model.Nt + time_iter) = u2(r_idx_1 * Ny + r_idx_2);
    }
}

void Acoustic::WriteData()
{
    json j_last_wavefield = Acoustic::EigenToJson(seis_data.last_wavefield);
    json j_seis_signal = Acoustic::EigenToJson(seis_data.seis_signal);
    json j_seis_signal_multi = Acoustic::EigenToJson(seis_data.seis_signal_multi);

    json output = {
        {"last_wavefield", j_last_wavefield},
        {"seis_signal", j_seis_signal},
        {"seis_signal_multi", j_seis_signal_multi}};

    std::ofstream o("./data/seis_data.json");
    o << std::setw(4) << output << std::endl;
}

json Acoustic::EigenToJson(const Eigen::VectorXd &vec)
{
    std::vector<double> v(vec.data(), vec.data() + vec.rows());
    json j_vec(v);
    return j_vec;
}
