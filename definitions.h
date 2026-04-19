#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define YPERIODIC
#define ZPERIODIC

#define LEFT_NEBB_VELOCITY
#define LEFT_U_VELOCITY 0.0
#define LEFT_V_VELOCITY 0.0
#define LEFT_W_VELOCITY 0.0

#define RIGHT_NEBB_VELOCITY
#define RIGHT_U_VELOCITY 0.0
#define RIGHT_V_VELOCITY 0.0
#define RIGHT_W_VELOCITY 0.0

// --------------------------
//     Initial conditions
// --------------------------
#define INI_TWOCOMPONENT_POISEUILLE
#define INI_TWOCOMPONENT_POISEUILLE_A   (0.25*width)
#define INI_TWOCOMPONENT_POISEUILLE_SF  2.0

#endif
