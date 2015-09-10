#include <fstream>
#include <sstream>
#include "SimulationParameters.hpp"
#include "WriteFiles.hpp"

/// *****************************
/// WRITE FLUID STATE TO VTK FILE
/// *****************************

// The fluid state is writen to a VTK file at each t_disk step.
// The following data is written:
// - density difference (density - 1)
// - x-component of velocity
// - y-component of velocity
// The following code is designed in such a way that the file can be read by ParaView.

void write_fluid_vtk(int time
  , double **density
  , double **velocity_x
  , double **velocity_y
  , double **force_x
  , double **force_y)
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

 for(int X = 0; X < Nx; ++X) {
   output_file << X + 0.5 << " ";
 }

 output_file << "\n";
 output_file << "Y_COORDINATES " << Ny - 2 << " float\n";

 for(int Y = 1; Y < Ny - 1; ++Y) {
   output_file << Y - 0.5 << " ";
 }

 output_file << "\n";
 output_file << "Z_COORDINATES " << 1 << " float\n";
 output_file << 0 << "\n";
 output_file << "POINT_DATA " << Nx * (Ny - 2) << "\n";

 /// Write density difference

 output_file << "SCALARS density_difference float 1\n";
 output_file << "LOOKUP_TABLE default\n";

 for(int Y = 1; Y < Ny - 1; ++Y) {
   for(int X = 0; X < Nx; ++X) {
     output_file << density[X][Y] - 1 << "\n";
   }
 }

 /// Write velocity

 output_file << "VECTORS velocity_vector float\n";

 for(int Y = 1; Y < Ny - 1; ++Y) {
   for(int X = 0; X < Nx; ++X) {
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

// The particle state (node positions) is writen to a VTK file at each t_disk step.
// The following code is designed in such a way that the file can be read by ParaView.

void write_particle_vtk(int time
  , particle_struct particle)
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

 for(int n = 0; n < particle_num_nodes; ++n) {
   output_file << particle.node[n].x << " " << particle.node[n].y << " 0\n";
 }

 /// Write lines between neighboring nodes

 output_file << "LINES " << particle_num_nodes << " " << 3 * particle_num_nodes << "\n";

 for(int n = 0; n < particle_num_nodes; ++n) {
   output_file << "2 " << n << " " << (n + 1) % particle_num_nodes << "\n";
 }

 /// Write vertices

 output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
 output_file << particle_num_nodes << " ";

 for(int n = 0; n < particle_num_nodes; ++n) {
   output_file << n << " ";
 }

 /// Close file

 output_file.close();

 return;
}
