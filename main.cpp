// This is an example 2D immersed boundary lattice Boltzmann method code.
// It uses the D2Q9 lattice including forcing term.
// Rigid bottom and top walls are parallel to the x-axis (channel).
// The flow is periodic along the x-axis.
// One particle is positioned in the flow.
// This particle can be rigid/deformable and stationary/moving.
//
// Last update 14-Mar-2012 by Timm Krueger.
// This code may be changed and distributed freely.
//
// The lattice velocities are defined according to the following scheme:
// index:   0  1  2  3  4  5  6  7  8
// ----------------------------------
// x:       0 +1 -1  0  0 +1 -1 +1 -1
// y:       0  0  0 +1 -1 +1 -1 -1 +1
//
// 8  3  5  ^y
//  \ | /   |   x
// 2- 0 -1   --->
//  / | \
// 6  4  7

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************
#include <cmath> // mathematical library
#include <cstdlib> // standard library
#include <fstream> // file streams
#include <iostream> // for the use of 'cout'
#include <sstream> // string streams
#include <vector> // vector containers
// square function; replaces SQ(x) by ((x) * (x)) in the code
#define SQ(x) ((x) * (x))

using namespace std; // permanently use the standard namespace

/// *********************
/// SIMULATION PARAMETERS
/// *********************
// These are the relevant simulation parameters.
// They can be changed by the user.
// If a bottom or top wall shall move in negative x-direction,
// a negative velocity has to be specified.
// Moving walls and gravity can be switched on simultaneously.
/// Simulation types
// Exactly one of the following options has to be defined.
// RIGID_CYLINDER
// - for simulation of, e.g., Karman vortex street
// - the cylinder is kept in space, it does not move
// DEFORMABLE_CYLINDER
// - for simulation of a deformable cylinder
// - the cylinder moves along with the flow
// DEFORMABLE_RBC
// - for simulation of a deformable red blood cell
// - the cell moves along with the flow

#define RIGID_CYLINDER
//#define DEFORMABLE_CYLINDER
//#define DEFORMABLE_RBC

const static double M_PI = 3.14159265;
/// Fluid/lattice properties
// number of lattice nodes along the x-axis (periodic)
const int Nx = 220;
// number of lattice nodes along the y-axis (including two wall nodes)
const int Ny = 62;
// relaxation time
const double tau = 0.6;
// number of time steps (running from 1 to t_num)
const int t_num = 50000;
// disk write time step (data will be written to the disk every t_disk step)
const int t_disk = 200;
// info time step (screen message will be printed every t_info step)
const int t_info = 1000;
// force density due to gravity (in positive x-direction)
const double gravity = 0.00000;
// velocity of the bottom wall (in positive x-direction)
const double wall_vel_bottom = 0.1;
// velocity of the top wall (in positive x-direction)
const double wall_vel_top = -wall_vel_bottom;
const bool can_move = false;

/// Particle properties
const int particle_num_nodes = 10; // number of surface nodes
const double particle_radius = 1; // radius
const double particle_stiffness = 2.0; // stiffness modulus
const double particle_bending = 0; // bending modulus
const double particle_center_x = Nx / 2; // center position (x-component)
const double particle_center_y = Ny / 4 - 1; // center position (y-component)
const double particle_torque = 0.0;  // particle initial torque;
double torque = 0.0;
double ang_vel = 0.0;

/// *****************
/// DECLARE VARIABLES
/// *****************
// relaxation frequency (inverse of relaxation time)
const double omega = 1. / tau;
double ***pop;
double ***pop_old; // LBM populations (old and new)
double **density; // fluid density
double **velocity_x; // fluid velocity (x-component)
double **velocity_y; // fluid velocity (y-component)
double **force_x; // fluid force (x-component)
double **force_y; // fluid force (y-component)
std::vector<double> pop_eq(9); // equilibrium populations
const std::vector<double> weight = {4.0 / 9.0,
    1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0}; // lattice weights

/// ******************
/// PARTICLE STRUCTURE
/// ******************
// The following code handles the object immersed in the flow.
// In the present implementation, only a single object can be put into the flow.
/// Structure for surface nodes
// Each node has a current x- and y-position and a reference x- and y-position.
struct Node
{
  /// Constructor
  Node() {
    x = 0;
    y = 0;
    x_ref = 0;
    y_ref = 0;
    vel_x = 0;
    vel_y = 0;
    force_x = 0;
    force_y = 0;
  }
  /// Elements
  double x; // current x-position
  double y; // current y-position
  double x_ref; // reference x-position
  double y_ref; // reference y-position
  double vel_x; // node velocity (x-component)
  double vel_y; // node velocity (y-component)
  double force_x; // node force (x-component)
  double force_y; // node force (y-component)
};

/// Structure for object (either cylinder or red blood cell)
struct Particle
{
  /// Constructor
  Particle() {
    num_nodes = particle_num_nodes;
    radius = particle_radius;
    stiffness = particle_stiffness;
    center.x = particle_center_x;
    center.y = particle_center_y;
    center.x_ref = particle_center_x;
    center.y_ref = particle_center_y;
    node = new Node[num_nodes];

    // The initial shape of the object is set in the following.
    // For a cylinder (rigid or deformable), the nodes define a circle.
    for (int n = 0; n < num_nodes; ++n) {
      auto theta = 2. * M_PI * (double) n / num_nodes;
      #if defined RIGID_CYLINDER || defined DEFORMABLE_CYLINDER
        node[n].x = center.x + radius * sin(theta);
        node[n].x_ref = center.x + radius * sin(theta);
        node[n].y = center.y + radius * cos(theta);
        node[n].y_ref = center.y + radius * cos(theta);
      #endif
      #ifdef DEFORMABLE_RBC
      #endif
    }  // n
  }
  /// Elements
  int num_nodes; // number of surface nodes
  double radius; // object radius
  double stiffness; // stiffness modulus
  Node center; // center node
  Node *node; // list of nodes
//  double torque;
};
/// *****************
/// DECLARE FUNCTIONS
/// *****************
// The following functions are used in the simulation code.
// allocate memory and initialize variables
void Initialize();
// perform LBM operations
void LatticeBoltzmannMethod(int time);
// compute fluid density and velocity from the populations
void ComputeMacroscopicProperties();
// compute the equilibrium populations from the fluid density and velocity
void ComputeEquilibrium(double den
  , double vel_x
  , double vel_y);
// compute the forces acting on the object nodes
void ComputeParticleForces(Particle particle
  , int time);
// spread node forces to fluid lattice
void SpreadForce(Particle particle);
// interpolate node velocities from fluid velocity
void InterpolateVelocity(Particle particle);
// update object center position
void UpdateParticlePosition(Particle particle);
// write the fluid state to the disk as VTK file
void WriteFluidVTK(int time);
// write the particle state to the disk as VTK file
void WriteParticleVTK(int time
  , Particle particle);
// write data to the disk (drag/lift, center position)
void WriteData(int time
  , Particle particle);

/// *************
/// MAIN FUNCTION
/// *************
// This is the main function, containing the simulation initialization and the
// simulation loop.
int main()
{
  /// ************
  /// PREPARATIONS
  /// ************
  Initialize(); // allocate memory and initialize variables
  Particle particle; // create immersed object
  /// Compute derived quantities
  // inner channel diameter
  const double D = Ny - 2;
  // lattice viscosity
  const double nu = (tau - 0.5) / 3;
  // expected maximum velocity for Poiseuille flow without immersed object
  const double umax = gravity / (2 * nu) * SQ(0.5 * D);
  // Reynolds number for Poiseuille flow without immersed object
  const double Re = D * umax / nu;

  /// Report derived parameters
  std::cout << "simulation parameters" << std::endl;
  std::cout << "=====================" << std::endl;
  std::cout << "D = " << D << std::endl;
  std::cout << "nu = " << nu << std::endl;
  std::cout << "umax = " << umax << std::endl;
  std::cout << "Re = " << Re << std::endl;
  std::cout << std::endl;
  /// ***************
  /// SIMULATION LOOP
  /// ***************
  // Overview of simulation algorithm:
  // 1) compute the node forces based on the object's deformation
  // 2) spread the node forces to the fluid lattice
  // 3) update the fluid state via LBM
  // 4) interpolate the fluid velocity to the object nodes
  // 5) update node positions and object's center
  // 6) if desired, write data to disk and report status
  std::cout << "starting simulation" << std::endl;
//  srand(1);
  // run over all times between 1 and t_num
  for (int t = 1; t <= t_num; ++t) {
    ComputeParticleForces(particle, t);
    SpreadForce(particle);
    LatticeBoltzmannMethod(t);
    InterpolateVelocity(particle);
    UpdateParticlePosition(particle);
    /// Write fluid and particle to VTK files
    // The data is only written each t_info time step.
    if (t % t_disk == 0) {
      WriteFluidVTK(t);
      WriteParticleVTK(t, particle);
      WriteData(t, particle);
    }
    if(t % t_info == 0) std::cout << t << " in " << t_num << std::endl;
  }  // t
  std::cout << "simulation complete" << std::endl;
  return 0;
} // end of main function

/// ****************************************
/// ALLOCATE MEMORY AND INITIALIZE VARIABLES
/// ****************************************
// The memory for lattice variables (populations, density, velocity, forces) is
// allocated. The variables are initialized.
void Initialize()
{
  /// Create folders, delete data file
  // Make sure that the VTK folders exist.
  // Old file data.dat is deleted, if existing.
  int ignore; // ignore return value of system calls
  ignore = system("mkdir -p vtk_fluid"); // create folder if not existing
  ignore = system("mkdir -p vtk_particle"); // create folder if not existing
  ignore = system("rm -f data.dat"); // delete file if existing

  /// Allocate memory for the fluid density, velocity, and force
  density = new double*[Nx];
  velocity_x = new double*[Nx];
  velocity_y = new double*[Nx];
  force_x = new double*[Nx];
  force_y = new double*[Nx];
  for (int X = 0; X < Nx; ++X) {
    density[X] = new double[Ny];
    velocity_x[X] = new double[Ny];
    velocity_y[X] = new double[Ny];
    force_x[X] = new double[Ny];
    force_y[X] = new double[Ny];
  }  // X
  /// Initialize the fluid density and velocity
  // Start with unit density and zero velocity.
  for (int X = 0; X < Nx; ++X) {
    for (int Y = 0; Y < Ny; ++Y) {
      density[X][Y] = 1;
      velocity_x[X][Y] = 0;
      velocity_y[X][Y] = 0;
      force_x[X][Y] = 0;
      force_y[X][Y] = 0;
    }  // Y
  }  // X
  /// Allocate memory for the populations
  pop = new double**[9];
  pop_old = new double**[9];
  for (int i = 0; i < 9; ++i) {
    pop[i] = new double*[Nx];
    pop_old[i] = new double*[Nx];
    for (int X = 0; X < Nx; ++X) {
      pop[i][X] = new double[Ny];
      pop_old[i][X] = new double[Ny];
      for (int Y = 0; Y < Ny; ++Y) {
        pop[i][X][Y] = 0;
        pop_old[i][X][Y] = 0;
      }  // Y
    }  // X
  }  // i
  /// Initialize the populations
  // Use the equilibrium populations corresponding to the initialized fluid
  // density and velocity.
  for (int X = 0; X < Nx; ++X) {
    for (int Y = 0; Y < Ny; ++Y) {
      ComputeEquilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
      for (int i = 0; i < 9; ++i) {
        pop_old[i][X][Y] = pop_eq[i];
        pop[i][X][Y] = pop_eq[i];
      }  // i
    }  // Y
  }  // X
  return;
}

/// *******************
/// COMPUTE EQUILIBRIUM
/// *******************
// This function computes the equilibrium populations from the fluid density and velocity.
// It computes the equilibrium only at a specific lattice node: Function has to be called at each lattice node.
// The standard quadratic euilibrium is used.
void ComputeEquilibrium(double den
  , double vel_x
  , double vel_y)
{
  const auto u_sqr = 1.5 * SQ(vel_x) + SQ(vel_y);
  pop_eq[0] = weight[0] * den * (1                                                     - u_sqr);
  pop_eq[1] = weight[1] * den * (1 + 3 * (  vel_x        ) + 4.5 * SQ(  vel_x        ) - u_sqr);
  pop_eq[2] = weight[2] * den * (1 + 3 * (- vel_x        ) + 4.5 * SQ(- vel_x        ) - u_sqr);
  pop_eq[3] = weight[3] * den * (1 + 3 * (          vel_y) + 4.5 * SQ(          vel_y) - u_sqr);
  pop_eq[4] = weight[4] * den * (1 + 3 * (        - vel_y) + 4.5 * SQ(        - vel_y) - u_sqr);
  pop_eq[5] = weight[5] * den * (1 + 3 * (  vel_x + vel_y) + 4.5 * SQ(  vel_x + vel_y) - u_sqr);
  pop_eq[6] = weight[6] * den * (1 + 3 * (- vel_x - vel_y) + 4.5 * SQ(- vel_x - vel_y) - u_sqr);
  pop_eq[7] = weight[7] * den * (1 + 3 * (  vel_x - vel_y) + 4.5 * SQ(  vel_x - vel_y) - u_sqr);
  pop_eq[8] = weight[8] * den * (1 + 3 * (- vel_x + vel_y) + 4.5 * SQ(- vel_x + vel_y) - u_sqr);

  return;
}

/// **********************
/// PERFORM LBM OPERATIONS
/// **********************

void LatticeBoltzmannMethod(int time)
{
  /// Swap populations
  // The present code used old and new populations which are swapped at the beginning of each time step.
  // This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
  // This way, the old populations are not overwritten during propagation.
  // The resulting code is easier to write and to debug.
  // The memory requirement for the populations is twice as large.
  double ***swap_temp = pop_old;
  pop_old = pop;
  pop = swap_temp;
  /// Lattice Boltzmann equation
  // The lattice Boltzmann equation is solved in the following.
  // The algorithm includes
  // - computation of the lattice force
  // - combined collision and propagation (faster than first collision and then propagation)
  for (int X = 0; X < Nx; ++X) {
    for (int Y = 1; Y < Ny - 1; ++Y) {
      /// Compute equilibrium
      // The equilibrium populations are computed.
      // Forces are coupled via Shan-Chen velocity shift.
      ComputeEquilibrium(density[X][Y]
        , (velocity_x[X][Y] + (force_x[X][Y] + gravity) * tau / density[X][Y])
        , (velocity_y[X][Y] + (force_y[X][Y]) * tau / density[X][Y]));
      /// Compute new populations
      // This is the lattice Boltzmann equation (combined collision and
      // propagation) including external forcing. Periodicity of the lattice in
      // x-direction is taken into account by the %-operator.
      pop[0][X]                [Y]     = pop_old[0][X][Y] * (1 - omega) + pop_eq[0] * omega;
      pop[1][(X + 1) % Nx]     [Y]     = pop_old[1][X][Y] * (1 - omega) + pop_eq[1] * omega;
      pop[2][(X - 1 + Nx) % Nx][Y]     = pop_old[2][X][Y] * (1 - omega) + pop_eq[2] * omega;
      pop[3][X]                [Y + 1] = pop_old[3][X][Y] * (1 - omega) + pop_eq[3] * omega;
      pop[4][X]                [Y - 1] = pop_old[4][X][Y] * (1 - omega) + pop_eq[4] * omega;
      pop[5][(X + 1) % Nx]     [Y + 1] = pop_old[5][X][Y] * (1 - omega) + pop_eq[5] * omega;
      pop[6][(X - 1 + Nx) % Nx][Y - 1] = pop_old[6][X][Y] * (1 - omega) + pop_eq[6] * omega;
      pop[7][(X + 1) % Nx]     [Y - 1] = pop_old[7][X][Y] * (1 - omega) + pop_eq[7] * omega;
      pop[8][(X - 1 + Nx) % Nx][Y + 1] = pop_old[8][X][Y] * (1 - omega) + pop_eq[8] * omega;
    }  // Y
  }  // X
  /// Bounce-back
  // Due to the presence of the rigid walls at y = 0 and y = Ny - 1, the
  // populations have to be bounced back. Ladd's momentum correction term is
  // included for moving walls (wall velocity parallel to x-axis). Periodicity
  // of the lattice in x-direction is taken into account via the %-operator.
  for (int X = 0; X < Nx; ++X) {
    /// Bottom wall (y = 0)
    pop[3][X][1] = pop[4][X]                [0];
    pop[5][X][1] = pop[6][(X - 1 + Nx) % Nx][0] + 6 * weight[6] * density[X][1] * wall_vel_bottom;
    pop[8][X][1] = pop[7][(X + 1) % Nx]     [0] - 6 * weight[7] * density[X][1] * wall_vel_bottom;
    /// Top wall (y = Ny - 1)
    pop[4][X][Ny - 2] = pop[3][X]                [Ny - 1];
    pop[6][X][Ny - 2] = pop[5][(X + 1) % Nx]     [Ny - 1] - 6 * weight[5] * density[X][Ny - 2] * wall_vel_top;
    pop[7][X][Ny - 2] = pop[8][(X - 1 + Nx) % Nx][Ny - 1] + 6 * weight[8] * density[X][Ny - 2] * wall_vel_top;
  }  // X
  ComputeMacroscopicProperties();
  return;
}

/// **********************************
/// COMPUTE FLUID DENSITY AND VELOCITY
/// **********************************
// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is not included here.
// It must be taken into account whenever the physical velocity is required.
void ComputeMacroscopicProperties()
{
  for (int X = 0; X < Nx; ++X) {
    for (int Y = 1; Y < Ny - 1; ++Y) {
      density[X][Y] = pop[0][X][Y] + pop[1][X][Y] + pop[2][X][Y] + pop[3][X][Y] + pop[4][X][Y] + pop[5][X][Y] + pop[6][X][Y] + pop[7][X][Y] + pop[8][X][Y];
      velocity_x[X][Y] = (pop[1][X][Y] - pop[2][X][Y] + pop[5][X][Y] - pop[6][X][Y] + pop[7][X][Y] - pop[8][X][Y]) / density[X][Y];
      velocity_y[X][Y] = (pop[3][X][Y] - pop[4][X][Y] + pop[5][X][Y] - pop[6][X][Y] - pop[7][X][Y] + pop[8][X][Y]) / density[X][Y];
    }
  }
  return;
}

/// ***********************
/// COMPUTE PARTICLE FORCES
/// ***********************
// The forces acting on the object nodes are computed.
// Depending on the simulation type (rigid/deformable cylinder or red blood cell),
// the force computation is different.
void ComputeParticleForces(Particle particle, int time)
{
  /// Reset forces
  // This way, the force from the previous time step is deleted.
  // This is necessary whenever forces are computed using '+='.
  for (int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].force_x = 0;
    particle.node[n].force_y = 0;
  }
  torque = 0.0;

  /// Compute strain forces
  #if defined DEFORMABLE_CYLINDER || defined DEFORMABLE_RBC
  #endif
  /// Compute bending forces
  #if defined DEFORMABLE_CYLINDER || defined DEFORMABLE_RBC
  #endif
  /// Compute rigid forces
  // Here, the node forces are proportional to the displacement with respect to the reference position.
  #ifdef RIGID_CYLINDER
    const double area = 2 * M_PI * particle.radius / particle.num_nodes; // area belonging to a node
    for (int n = 0; n < particle.num_nodes; ++n) {
      auto restore_force_x = -particle.stiffness * (particle.node[n].x - particle.node[n].x_ref) * area;
      auto restore_force_y = -particle.stiffness * (particle.node[n].y - particle.node[n].y_ref) * area;
      particle.node[n].force_x += restore_force_x;
      particle.node[n].force_y += restore_force_y;
      const auto x_d = particle.node[n].x - particle.center.x;
      const auto y_d = particle.node[n].y - particle.center.y;
      torque += x_d * -restore_force_y - y_d * -restore_force_x;
    }
//    std::cout << torque << std::endl;
  #endif
  return;
}

/// *************
/// SPREAD FORCES
/// *************
// The node forces are spread to the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the
// present code.
void SpreadForce(Particle particle)
{
  /// Reset forces
  // This is necessary since '+=' is used afterwards.
  for (int X = 0; X < Nx; ++X) {
    for (int Y = 1; Y < Ny - 1; ++Y) {
      force_x[X][Y] = 0;
      force_y[X][Y] = 0;
    }
  }
  /// Spread forces
  // Run over all object nodes.
  for (int n = 0; n < particle.num_nodes; ++n) {
    // Identify the lowest fluid lattice node in interpolation range.
    // 'Lowest' means: its x- and y-values are the smallest.
    // The other fluid nodes in range have coordinates
    // (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).
    int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
    int y_int = (int) (particle.node[n].y + 0.5 );

    // Run over all neighboring fluid nodes.
    // In the case of the two-point interpolation, it is 2x2 fluid nodes.
    for (int X = x_int; X <= x_int + 1; ++X) {
      for (int Y = y_int; Y <= y_int + 1; ++Y) {
        // Compute distance between object node and fluid lattice node.
        const double dist_x = particle.node[n].x - 0.5 - X;
        const double dist_y = particle.node[n].y + 0.5 - Y;
        // Compute interpolation weights for x- and y-direction based on the distance.
        const double weight_x = 1 - abs(dist_x);
        const double weight_y = 1 - abs(dist_y);
        const auto weights = weight_x * weight_y;
        // Compute lattice force.
        force_x[(X + Nx) % Nx][Y] += particle.node[n].force_x * weights;
        force_y[(X + Nx) % Nx][Y] += particle.node[n].force_y * weights;
      }
    }
  }
  return;
}

/// **********************
/// INTERPOLATE VELOCITIES
/// **********************
// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.
void InterpolateVelocity(Particle particle)
{
  // Run over all object nodes.
  for (int n = 0; n < particle.num_nodes; ++n) {
    // Reset node velocity first since '+=' is used.
    particle.node[n].vel_x = 0;
    particle.node[n].vel_y = 0;
    // Identify the lowest fluid lattice node in interpolation range (see spreading).
    int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
    int y_int = (int) (particle.node[n].y + 0.5 );
    // Run over all neighboring fluid nodes.
    // In the case of the two-point interpolation, it is 2x2 fluid nodes.
    for (int X = x_int; X <= x_int + 1; ++X) {
      for (int Y = y_int; Y <= y_int + 1; ++Y) {
        // Compute distance between object node and fluid lattice node.
        const double dist_x = particle.node[n].x - 0.5 - X;
        const double dist_y = particle.node[n].y + 0.5 - Y;
        // Compute interpolation weights for x- and y-direction based on the distance.
        const double weight_x = 1 - abs(dist_x);
        const double weight_y = 1 - abs(dist_y);
        const auto weights = weight_x * weight_y;
        // Compute node velocities.
        particle.node[n].vel_x += ((velocity_x[(X + Nx) % Nx][Y] + 0.5 * (force_x[(X + Nx) % Nx][Y] + gravity) / density[(X + Nx) % Nx][Y]) * weights);
        particle.node[n].vel_y += ((velocity_y[(X + Nx) % Nx][Y] + 0.5 * (force_y[(X + Nx) % Nx][Y]) / density[(X + Nx) % Nx][Y]) * weights);
      }
    }
  }
  return;
}

/// ************************
/// UPDATE PARTICLE POSITION
/// ************************
// The position of the particle nodes are updated according to their velocity.
// The center position is updated as well.
// The new node position is its old position plus its current velocity (Euler
// integration). The center position is the arithmetic mean of all node
// positions. Periodicity is taken into account: If the particle center leaves
// the system domain (x < 0 or x >= Nx), it reenters from the other side.
void UpdateParticlePosition(Particle particle)
{
  /// Reset center position
  particle.center.x = 0;
  particle.center.y = 0;
  /// Update node and center positions
  for (int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].x += particle.node[n].vel_x;
    particle.node[n].y += particle.node[n].vel_y;
    particle.center.x += particle.node[n].x / particle.num_nodes;
    particle.center.y += particle.node[n].y / particle.num_nodes;
  }
  if (can_move) {
    particle.center.x_ref = particle.center.x;
    particle.center.y_ref = particle.center.y;
    for (int n = 0; n < particle.num_nodes; ++n) {
      particle.node[n].x_ref = particle.center.x_ref + particle.radius *
          sin(2. *  M_PI * (double) n / particle.num_nodes);
      particle.node[n].y_ref = particle.center.y_ref + particle.radius *
          cos(2. *  M_PI * (double) n / particle.num_nodes);
    }  // n
//    for (int n = 0; n < particle.num_nodes; ++n) {
//      particle.node[n].x += particle.node[n].vel_x;
//      particle.node[n].y += particle.node[n].vel_y;
//    }
  }
  ang_vel += torque / particle_radius / particle_radius;
  std::cout << ang_vel << std::endl;
//  const auto angle = fabs(ang_vel);
//  const auto dir = ang_vel / fabs(ang_vel);
  for (int n = 0; n < particle.num_nodes; ++n) {
    const auto old_x_ref = particle.node[n].x_ref - particle.center.x_ref;
    const auto old_y_ref = particle.node[n].y_ref - particle.center.y_ref;
    const auto old_x = particle.node[n].x - particle.center.x;
    const auto old_y = particle.node[n].y - particle.center.y;
    const auto mov_x_ref = cos(ang_vel) * old_x_ref - sin(ang_vel) * old_y_ref;
    const auto mov_y_ref = sin(ang_vel) * old_x_ref + cos(ang_vel) * old_y_ref;
    const auto mov_x = cos(ang_vel) * old_x - sin(ang_vel) * old_y;
    const auto mov_y = sin(ang_vel) * old_x + cos(ang_vel) * old_y;
//    std::cout << mov_x << std::endl;
    particle.node[n].x_ref = particle.center.x_ref + mov_x_ref;
    particle.node[n].y_ref = particle.center.y_ref + mov_y_ref;
    particle.node[n].x = particle.center.x + mov_x;
    particle.node[n].y = particle.center.y + mov_y;
  }  // n

  /// Check for periodicity along the x-axis
  if (particle.center.x < 0) {
    particle.center.x += Nx;
    for (int n = 0; n < particle.num_nodes; ++n) particle.node[n].x += Nx;
  }
  else if (particle.center.x >= Nx) {
    particle.center.x -= Nx;
    for (int n = 0; n < particle.num_nodes; ++n) particle.node[n].x -= Nx;
  }
  if (can_move) {
    if (particle.center.x_ref < 0) {
      particle.center.x_ref += Nx;
      for (int n = 0; n < particle.num_nodes; ++n) particle.node[n].x_ref += Nx;
    }
    else if (particle.center.x_ref >= Nx) {
      particle.center.x_ref -= Nx;
      for (int n = 0; n < particle.num_nodes; ++n) particle.node[n].x_ref -= Nx;
    }
  }
  return;
}

/// *****************************
/// WRITE FLUID STATE TO VTK FILE
/// *****************************
// The fluid state is writen to a VTK file at each t_disk step.
// The following data is written:
// - density difference (density - 1)
// - x-component of velocity
// - y-component of velocity
// The following code is designed in such a way that the file can be read by
// ParaView.
void WriteFluidVTK(int time)
{
  /// Create filename
  std::stringstream output_filename;
  output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
  std::ofstream output_file;

  /// Open file
  output_file.open(output_filename.str().c_str());
  /// Write VTK header
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "fluid_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET RECTILINEAR_GRID\n";
  output_file << "DIMENSIONS " << Nx << " " << Ny - 2 << " 1" << "\n";
  output_file << "X_COORDINATES " << Nx << " float\n";

  for (int X = 0; X < Nx; ++X) output_file << X + 0.5 << " ";
  output_file << "\n";
  output_file << "Y_COORDINATES " << Ny - 2 << " float\n";
  for (int Y = 1; Y < Ny - 1; ++Y) output_file << Y - 0.5 << " ";
  output_file << "\n";
  output_file << "Z_COORDINATES " << 1 << " float\n";
  output_file << 0 << "\n";
  output_file << "POINT_DATA " << Nx * (Ny - 2) << "\n";

  /// Write density difference
  output_file << "SCALARS density_difference float 1\n";
  output_file << "LOOKUP_TABLE default\n";

  for (int Y = 1; Y < Ny - 1; ++Y) {
    for (int X = 0; X < Nx; ++X) output_file << density[X][Y] - 1 << "\n";
  }  // Y

  /// Write velocity
  output_file << "VECTORS velocity_vector float\n";
  for (int Y = 1; Y < Ny - 1; ++Y) {
    for (int X = 0; X < Nx; ++X) {
      output_file << velocity_x[X][Y] + 0.5 * (force_x[X][Y] + gravity) / density[X][Y] << " " << velocity_y[X][Y] + 0.5 * (force_y[X][Y]) / density[X][Y] << " 0\n";
    }
  }
  /// Close file
  output_file.close();
  return;
}

/// ********************************
/// WRITE PARTICLE STATE TO VTK FILE
/// ********************************
// The particle state (node positions) is writen to a VTK file at each t_disk
// step. The following code is designed in such a way that the file can be read
// by ParaView.

void WriteParticleVTK(int time, Particle particle)
{
  /// Create filename
  std::stringstream output_filename;
  output_filename << "vtk_particle/particle_t" << time << ".vtk";
  std::ofstream output_file;

  /// Open file
  output_file.open(output_filename.str().c_str());

  /// Write VTK header
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << particle_num_nodes << " float\n";
  for (int n = 0; n < particle_num_nodes; ++n) {
    output_file << particle.node[n].x << " " << particle.node[n].y << " 0\n";
  }

  /// Write lines between neighboring nodes
  output_file << "LINES " << particle_num_nodes << " " << 3 * particle_num_nodes << "\n";
  for (int n = 0; n < particle_num_nodes; ++n) {
    output_file << "2 " << n << " " << (n + 1) % particle_num_nodes << "\n";
  }

  /// Write vertices
  output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
  output_file << particle_num_nodes << " ";
  for (int n = 0; n < particle_num_nodes; ++n) output_file << n << " ";

  /// Close file
  output_file.close();
  return;
}

/// ************************
/// WRITE DATA TO ASCII FILE
/// ************************
// The following quantities are written to the disk at each t_disk step:
// - drag and lift forces (x- and y-components of the force)
// - object center position (x- and y-components)
// - object center velocity (x- and y-components)
// The data file is readable by gnuplot
void WriteData(int time, Particle particle)
{
  /// Create filename
  std::string output_filename("data.dat");
  std::ofstream output_file;

  /// Open file
  output_file.open(output_filename.c_str(), std::fstream::app);

  /// Compute quantities
  double force_tot_x = 0;
  double force_tot_y = 0;
  double vel_center_x = 0;
  double vel_center_y = 0;

  for (int i = 0; i < particle.num_nodes; ++i) {
    force_tot_x += particle.node[i].force_x;
    force_tot_y += particle.node[i].force_y;
    vel_center_x += particle.node[i].vel_x;
    vel_center_y += particle.node[i].vel_y;
  }

  /// Write data
  output_file << time << " "; // time step
  output_file << force_tot_x << " "; // drag force
  output_file << force_tot_y << " "; // lift force
  output_file << particle.center.x << " "; // center position (x-component)
  output_file << particle.center.y << " "; // center position (y-component)
  output_file << particle.center.y_ref << " "; // center position (y-component)
  output_file << vel_center_x << " "; // center velocity (x-component)
  output_file << vel_center_y << "\n"; // center velocity (y-component)

  /// Close file
  output_file.close();
  return;
}
