#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define YPERIODIC
#define ZPERIODIC

#define LEFT_NEBB_VELOCITY
#define LEFT_U_VELOCITY 0.0
#define LEFT_V_VELOCITY 1e-2
#define LEFT_W_VELOCITY 0.0

#define RIGHT_NEBB_VELOCITY
#define RIGHT_U_VELOCITY 0.0
#define RIGHT_V_VELOCITY 0.0
#define RIGHT_W_VELOCITY 0.0

// --------------------------
//     Initial conditions
// --------------------------
#define INI_TWOCOMPONENT_COUETTE
#define INI_TWOCOMPONENT_COUETTE_V_LEFT 1e-2
#define INI_TWOCOMPONENT_COUETTE_SF     1.0

#endif
