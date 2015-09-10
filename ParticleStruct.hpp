#ifndef PARTICLE_STRUCT_HPP_
#define PARTICLE_STRUCT_HPP_
#include <cmath>

/// Particle properties

const int particle_num_nodes = 60; // number of surface nodes
const double particle_radius = 8; // radius
const double particle_stiffness = 1; // stiffness modulus
const double particle_bending = 0; // bending modulus
const double particle_center_x = 30; // center position (x-component)
const double particle_center_y = 15; // center position (y-component)

/// ******************
/// PARTICLE STRUCTURE
/// ******************

// The following code handles the object immersed in the flow.
// In the present implementation, only a single object can be put into the flow.

/// Structure for surface nodes
// Each node has a current x- and y-position and a reference x- and y-position.

struct node_struct {

 /// Constructor

 node_struct() {
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

struct particle_struct {

 /// Constructor

 particle_struct() {
   num_nodes = particle_num_nodes;
   radius = particle_radius;
   stiffness = particle_stiffness;
   u_x = 0.0;
   u_y = 0.0;
   center.x = particle_center_x;
   center.y = particle_center_y;
   center.x_ref = particle_center_x;
   center.y_ref = particle_center_y;
   node = new node_struct[num_nodes];

   // The initial shape of the object is set in the following.
   // For a cylinder (rigid or deformable), the nodes define a circle.
   // For a red blood cell, the y-position has to be changed in order to describe a red blood cell.
   // Initially, the current node positions and reference node positions are identical.
   // During the simulation, only the current positions are updated,
   // the reference node positions are fixed.

   for(int n = 0; n < num_nodes; ++n) {
     node[n].x = center.x + radius * sin(2. * M_PI * (double) n / num_nodes);
     node[n].x_ref = center.x + radius * sin(2. * M_PI * (double) n / num_nodes);
     node[n].y = center.y + radius * cos(2. * M_PI * (double) n / num_nodes);
     node[n].y_ref = center.y + radius * cos(2. * M_PI * (double) n / num_nodes);
   }
 }

 /// Elements

 int num_nodes; // number of surface nodes
 double radius; // object radius
 double stiffness; // stiffness modulus
 double u_x;  // x vel of particle
 double u_y;  // y vel of particle
 node_struct center; // center node
 node_struct *node; // list of nodes
};

#endif // PARTICLE_STRUCT_HPP_
