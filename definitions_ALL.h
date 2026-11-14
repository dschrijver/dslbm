#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//         Dimension
// ---------------------------
#define DIM_2D
#define DIM_3D

// ---------------------------
//          Stencil
// ---------------------------
#define D2Q9
#define D3Q27

// ---------------------------
//     Boundary conditions
// ---------------------------
#define XPERIODIC
#define YPERIODIC
#define ZPERIODIC

// --------------------------
//     Initial conditions
// --------------------------
#define INI_FLUID_1 // Stationary fluid, constant density, random orientation

#define INI_ACTIVE_1 // Random orientation

#endif
