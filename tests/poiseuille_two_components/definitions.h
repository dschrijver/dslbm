#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define YPERIODIC
#define ZPERIODIC

#define LEFT_HWBB_NOSLIP
#define RIGHT_HWBB_NOSLIP

// --------------------------
//     Initial conditions
// --------------------------
#define INI_TWOCOMPONENT_POISEUILLE
#define INI_TWOCOMPONENT_POISEUILLE_A   (0.25*width)
#define INI_TWOCOMPONENT_POISEUILLE_SF  2.0

#endif
