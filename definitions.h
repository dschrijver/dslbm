#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// #define NO_NEMATICS // Turns off all nematics related machinery.

// ---------------------------
//         Dimension
// ---------------------------
#define DIM_2D
// #define DIM_3D

// ---------------------------
//          Stencil
// ---------------------------
#define D2Q9
// #define D3Q27

// ---------------------------
//     Boundary conditions
// ---------------------------
#define XPERIODIC
#define YPERIODIC
#define ZPERIODIC

// #define LEFT_NOSLIP
// #define RIGHT_NOSLIP
// #define BOTTOM_NOSLIP
// #define TOP_NOSLIP
// #define BACK_NOSLIP
// #define FRONT_NOSLIP

// --------------------------
//     Initial conditions
// --------------------------
// Fluid initial conditions
#define INI_FLUID_1 // Stationary fluid, constant density

// Activity initial conditions
#define INI_ACTIVE_1 // Random orientation

#endif
