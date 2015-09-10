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
// 8 3 5  ^y
//  \|/   |   x
// 2-0-1   --->
//  /|\
// 6 4 7

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************

#include <vector> // vector containers
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library
#include "ParticleStruct.hpp"
#include "SimulationParameters.hpp"
#include "WriteFiles.hpp"
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

using namespace std; // permanently use the standard namespace

/// *********************
/// SIMULATION PARAMETERS
/// *********************

// These are the relevant simulation parameters.
// They can be changed by the user.
// If a bottom or top wall shall move in negative x-direction, a negative velocity has to be specified.
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

/// *****************
/// DECLARE VARIABLES
/// *****************

// The following code should not be modified when it is first used.

const double omega = 1. / tau; // relaxation frequency (inverse of relaxation time)
double ***pop, ***pop_old; // LBM populations (old and new)
double **density; // fluid density
double **velocity_x; // fluid velocity (x-component)
double **velocity_y; // fluid velocity (y-component)
double **force_x; // fluid force (x-component)
double **force_y; // fluid force (y-component)
double pop_eq[9]; // equilibrium populations
const double weight[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}; // lattice weights

/// *****************
/// DECLARE FUNCTIONS
/// *****************

// The following functions are used in the simulation code.

void initialize(); // allocate memory and initialize variables
void LBM(int); // perform LBM operations
void momenta(); // compute fluid density and velocity from the populations
void equilibrium(double, double, double); // compute the equilibrium populations from the fluid density and velocity
void compute_particle_forces(particle_struct, int); // compute the forces acting on the object nodes
void spread(particle_struct); // spread node forces to fluid lattice
void interpolate(particle_struct); // interpolate node velocities from fluid velocity
void update_particle_position(particle_struct); // update object center position


/// *************
/// MAIN FUNCTION
/// *************

// This is the main function, containing the simulation initialization and the simulation loop.

int main() {

  /// ************
  /// PREPARATIONS
  /// ************

  initialize(); // allocate memory and initialize variables
  particle_struct particle; // create immersed object

  /// Compute derived quantities

  const double D = Ny - 2; // inner channel diameter
  const double nu = (tau - 0.5) / 3; // lattice viscosity
  const double umax = gravity / (2 * nu) * SQ(0.5 * D); // expected maximum velocity for Poiseuille flow without immersed object
  const double Re = D * umax / nu; // Reynolds number for Poiseuille flow without immersed object

  /// Report derived parameters

  cout << "simulation parameters" << endl;
  cout << "=====================" << endl;
  cout << "D = " << D << endl;
  cout << "nu = " << nu << endl;
  cout << "umax = " << umax << endl;
  cout << "Re = " << Re << endl;
  cout << endl;

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

  cout << "starting simulation" << endl;

  srand(1);

  for(int t = 1; t <= t_num; ++t) { // run over all times between 1 and t_num

    compute_particle_forces(particle, t); // compute particle forces
    spread(particle); // spread forces from the Lagrangian to the Eulerian mesh
    LBM(t); // perform collision, propagation, and bounce-back
    interpolate(particle); // interpolate velocity
    update_particle_position(particle); // update particle position

    /// Write fluid and particle to VTK files
    // The data is only written each t_info time step.

    if(t % t_disk == 0) {
      write_fluid_vtk(t, density, velocity_x, velocity_y, force_x, force_y);
      write_particle_vtk(t, particle);
    }

   /// Report end of time step

    if(t % t_info == 0) {
      cout << "completed time step " << t << " in [1, " << t_num << "]" << endl;
    }
  }  // t

  /// Report successful end of simulation

  cout << "simulation complete" << endl;

  return 0;
} // end of main function

/// ****************************************
/// ALLOCATE MEMORY AND INITIALIZE VARIABLES
/// ****************************************

// The memory for lattice variables (populations, density, velocity, forces) is allocated.
// The variables are initialized.

void initialize() {

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

  for(int X = 0; X < Nx; ++X) {
    density[X] = new double[Ny];
    velocity_x[X] = new double[Ny];
    velocity_y[X] = new double[Ny];
    force_x[X] = new double[Ny];
    force_y[X] = new double[Ny];
  }  // X

  /// Initialize the fluid density and velocity
  // Start with unit density and zero velocity.

  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
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

  for (int c_i = 0; c_i < 9; ++c_i) {
    pop[c_i] = new double*[Nx];
    pop_old[c_i] = new double*[Nx];

    for (int X = 0; X < Nx; ++X) {
      pop[c_i][X] = new double[Ny];
      pop_old[c_i][X] = new double[Ny];

      for (int Y = 0; Y < Ny; ++Y) {
        pop[c_i][X][Y] = 0;
        pop_old[c_i][X][Y] = 0;
      }  // Y
    }  // X
  }  // c_i

  /// Initialize the populations
  // Use the equilibrium populations corresponding to the initialized fluid density and velocity.

  for (int X = 0; X < Nx; ++X) {
    for (int Y = 0; Y < Ny; ++Y) {
      equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);

      for (int c_i = 0; c_i < 9; ++c_i) {
        pop_old[c_i][X][Y] = pop_eq[c_i];
        pop[c_i][X][Y] = pop_eq[c_i];
      }  // c_i
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
// reminder: SQ(x) = x * x

void equilibrium(double den, double vel_x, double vel_y) {
  pop_eq[0] = weight[0] * den * (1                                                     - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[1] = weight[1] * den * (1 + 3 * (  vel_x        ) + 4.5 * SQ(  vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[2] = weight[2] * den * (1 + 3 * (- vel_x        ) + 4.5 * SQ(- vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[3] = weight[3] * den * (1 + 3 * (          vel_y) + 4.5 * SQ(          vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[4] = weight[4] * den * (1 + 3 * (        - vel_y) + 4.5 * SQ(        - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[5] = weight[5] * den * (1 + 3 * (  vel_x + vel_y) + 4.5 * SQ(  vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[6] = weight[6] * den * (1 + 3 * (- vel_x - vel_y) + 4.5 * SQ(- vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[7] = weight[7] * den * (1 + 3 * (  vel_x - vel_y) + 4.5 * SQ(  vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[8] = weight[8] * den * (1 + 3 * (- vel_x + vel_y) + 4.5 * SQ(- vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));

  return;
}

/// **********************
/// PERFORM LBM OPERATIONS
/// **********************

void LBM(int time) {

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

  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {

      /// Compute equilibrium
      // The equilibrium populations are computed.
      // Forces are coupled via Shan-Chen velocity shift.

      equilibrium(density[X][Y], (velocity_x[X][Y] + (force_x[X][Y] + gravity) * tau / density[X][Y]), (velocity_y[X][Y] + (force_y[X][Y]) * tau / density[X][Y]));

      /// Compute new populations
      // This is the lattice Boltzmann equation (combined collision and propagation) including external forcing.
      // Periodicity of the lattice in x-direction is taken into account by the %-operator.

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
  // Due to the presence of the rigid walls at y = 0 and y = Ny - 1, the populations have to be bounced back.
  // Ladd's momentum correction term is included for moving walls (wall velocity parallel to x-axis).
  // Periodicity of the lattice in x-direction is taken into account via the %-operator.

  for(int X = 0; X < Nx; ++X) {

    /// Bottom wall (y = 0)

    pop[3][X][1] = pop[4][X]                [0];
    pop[5][X][1] = pop[6][(X - 1 + Nx) % Nx][0] + 6 * weight[6] * density[X][1] * wall_vel_bottom;
    pop[8][X][1] = pop[7][(X + 1) % Nx]     [0] - 6 * weight[7] * density[X][1] * wall_vel_bottom;

    /// Top wall (y = Ny - 1)

    pop[4][X][Ny - 2] = pop[3][X]                [Ny - 1];
    pop[6][X][Ny - 2] = pop[5][(X + 1) % Nx]     [Ny - 1] - 6 * weight[5] * density[X][Ny - 2] * wall_vel_top;
    pop[7][X][Ny - 2] = pop[8][(X - 1 + Nx) % Nx][Ny - 1] + 6 * weight[8] * density[X][Ny - 2] * wall_vel_top;
  }  // X

  /// Compute fluid density and velocity
  // The fluid density and velocity are obtained from the populations.

  momenta();

  return;
}

/// **********************************
/// COMPUTE FLUID DENSITY AND VELOCITY
/// **********************************

// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is not included here.
// It must be taken into account whenever the physical velocity is required.

void momenta() {
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {
      density[X][Y] = pop[0][X][Y] + pop[1][X][Y] + pop[2][X][Y] + pop[3][X][Y] + pop[4][X][Y] + pop[5][X][Y] + pop[6][X][Y] + pop[7][X][Y] + pop[8][X][Y];
      velocity_x[X][Y] = (pop[1][X][Y] - pop[2][X][Y] + pop[5][X][Y] - pop[6][X][Y] + pop[7][X][Y] - pop[8][X][Y]) / density[X][Y];
      velocity_y[X][Y] = (pop[3][X][Y] - pop[4][X][Y] + pop[5][X][Y] - pop[6][X][Y] - pop[7][X][Y] + pop[8][X][Y]) / density[X][Y];
    }  // Y
  }  // X

  return;
}

/// ***********************
/// COMPUTE PARTICLE FORCES
/// ***********************

// The forces acting on the object nodes are computed.
// Depending on the simulation type (rigid/deformable cylinder or red blood cell),
// the force computation is different.

void compute_particle_forces(particle_struct particle, int time) {

  /// Reset forces
  // This way, the force from the previous time step is deleted.
  // This is necessary whenever forces are computed using '+='.

  for(int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].force_x = 0;
    particle.node[n].force_y = 0;
  }  // n

  /// Compute rigid forces
  // Here, the node forces are proportional to the displacement with respect to the reference position.

  const double area = 2 * M_PI * particle.radius / particle.num_nodes; // area belonging to a node

  for(int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].force_x += -particle.stiffness * (particle.node[n].x - particle.node[n].x_ref) * area;
    particle.node[n].force_y += -particle.stiffness * (particle.node[n].y - particle.node[n].y_ref) * area;
  }  // n

  return;
}

/// *************
/// SPREAD FORCES
/// *************

// The node forces are spread to the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.

void spread(particle_struct particle) {

  /// Reset forces
  // This is necessary since '+=' is used afterwards.

  for(int X = 0; X < Nx; ++X) {
   for(int Y = 1; Y < Ny - 1; ++Y) {
     force_x[X][Y] = 0;
     force_y[X][Y] = 0;
   }
  }

  /// Spread forces
  // Run over all object nodes.

  for(int n = 0; n < particle.num_nodes; ++n) {

   // Identify the lowest fluid lattice node in interpolation range.
   // 'Lowest' means: its x- and y-values are the smallest.
   // The other fluid nodes in range have coordinates
   // (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).

   int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
   int y_int = (int) (particle.node[n].y + 0.5 );

   // Run over all neighboring fluid nodes.
   // In the case of the two-point interpolation, it is 2x2 fluid nodes.

   for(int X = x_int; X <= x_int + 1; ++X) {
      for(int Y = y_int; Y <= y_int + 1; ++Y) {

        // Compute distance between object node and fluid lattice node.

        const double dist_x = particle.node[n].x - 0.5 - X;
        const double dist_y = particle.node[n].y + 0.5 - Y;

        // Compute interpolation weights for x- and y-direction based on the distance.

        const double weight_x = 1 - abs(dist_x);
        const double weight_y = 1 - abs(dist_y);

        // Compute lattice force.

        force_x[(X + Nx) % Nx][Y] += (particle.node[n].force_x * weight_x * weight_y);
        force_y[(X + Nx) % Nx][Y] += (particle.node[n].force_y * weight_x * weight_y);
      }  // Y
    }  // X
  }  // n

  return;
}

/// **********************
/// INTERPOLATE VELOCITIES
/// **********************

// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.

void interpolate(particle_struct particle) {

  // Run over all object nodes.

  for(int n = 0; n < particle.num_nodes; ++n) {

   // Reset node velocity first since '+=' is used.

   particle.node[n].vel_x = 0;
   particle.node[n].vel_y = 0;

   // Identify the lowest fluid lattice node in interpolation range (see spreading).

   int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
   int y_int = (int) (particle.node[n].y + 0.5 );

   // Run over all neighboring fluid nodes.
   // In the case of the two-point interpolation, it is 2x2 fluid nodes.

   for(int X = x_int; X <= x_int + 1; ++X) {
      for(int Y = y_int; Y <= y_int + 1; ++Y) {

        // Compute distance between object node and fluid lattice node.

        const double dist_x = particle.node[n].x - 0.5 - X;
        const double dist_y = particle.node[n].y + 0.5 - Y;

        // Compute interpolation weights for x- and y-direction based on the distance.

        const double weight_x = 1 - abs(dist_x);
        const double weight_y = 1 - abs(dist_y);

        // Compute node velocities.

        particle.node[n].vel_x += ((velocity_x[(X + Nx) % Nx][Y] + 0.5 * (force_x[(X + Nx) % Nx][Y] + gravity) / density[(X + Nx) % Nx][Y]) * weight_x * weight_y);
        particle.node[n].vel_y += ((velocity_y[(X + Nx) % Nx][Y] + 0.5 * (force_y[(X + Nx) % Nx][Y]) / density[(X + Nx) % Nx][Y]) * weight_x * weight_y);
      }  // Y
    }  // X
  }  // n

  return;
}

/// ************************
/// UPDATE PARTICLE POSITION
/// ************************

// The position of the particle nodes are updated according to their velocity.
// The center position is updated as well.
// The new node position is its old position plus its current velocity (Euler integration).
// The center position is the arithmetic mean of all node positions.
// Periodicity is taken into account:
// If the particle center leaves the system domain (x < 0 or x >= Nx), it reenters from the other side.

void update_particle_position(particle_struct particle) {

  /// Reset center position

  particle.center.x = 0;
  particle.center.y = 0;

  /// Update node and center positions

  for(int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].x += particle.node[n].vel_x;
    particle.node[n].y += particle.node[n].vel_y;
    particle.center.x += particle.node[n].x / particle.num_nodes;
    particle.center.y += particle.node[n].y / particle.num_nodes;
  }

  /// Check for periodicity along the x-axis

  if(particle.center.x < 0) {
    particle.center.x += Nx;
    for(int n = 0; n < particle.num_nodes; ++n) particle.node[n].x += Nx;
  }
  else if(particle.center.x >= Nx) {
    particle.center.x -= Nx;
    for (int n = 0; n < particle.num_nodes; ++n) particle.node[n].x -= Nx;
  }

  return;
}
