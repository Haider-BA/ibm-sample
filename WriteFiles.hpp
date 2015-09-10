#ifndef WRITE_FILES_HPP_
#define WRITE_FILES_HPP_
#include "ParticleStruct.hpp"

void write_fluid_vtk(int time
  , double **density
  , double **velocity_x
  , double **velocity_y
  , double **force_x
  , double **force_y); // write the fluid state to the disk as VTK file

void write_particle_vtk(int time
  , particle_struct particle); // write the particle state to the disk as VTK file
#endif  // WRITE_FILES_HPP_
