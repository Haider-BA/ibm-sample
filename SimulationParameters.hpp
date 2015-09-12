#ifndef SIMULATION_PARAMETERS_HPP_
#define SIMULATION_PARAMETERS_HPP_

/// Fluid/lattice properties

const int Nx = 220; // number of lattice nodes along the x-axis (periodic)
const int Ny = 62; // number of lattice nodes along the y-axis (including two wall nodes)
const double tau = 1.0; // relaxation time
const int t_num = 50000; // number of time steps (running from 1 to t_num)
const int t_disk = 200; // disk write time step (data will be written to the disk every t_disk step)
const int t_info = 1000; // info time step (screen message will be printed every t_info step)
const double gravity = 0.0; // force density due to gravity (in positive x-direction)
const double wall_vel_bottom = 0.2; // velocity of the bottom wall (in positive x-direction)
const double wall_vel_top = -0.2; // velocity of the top wall (in positive x-direction)
#endif // SIMULATION_PARAMETERS_HPP_
