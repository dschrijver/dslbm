#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define XPERIODIC
#define YPERIODIC
#define ZPERIODIC

// --- Wetting ---
#define THETA_C 60.0

// --- Half-way Bounce-Back (HWBB) ---
#define LEFT_HWBB_NOSLIP
#define RIGHT_HWBB_NOSLIP
#define BOTTOM_HWBB_NOSLIP
#define TOP_HWBB_NOSLIP
#define BACK_HWBB_NOSLIP
#define FRONT_HWBB_NOSLIP

// --- Non-Equilibrium Bounce-Back ---
#define LEFT_NEBB_VELOCITY
#define LEFT_U_VELOCITY 0.0
#define LEFT_V_VELOCITY 0.0
#define LEFT_W_VELOCITY 0.0

#define RIGHT_NEBB_VELOCITY
#define RIGHT_U_VELOCITY 0.0
#define RIGHT_V_VELOCITY 0.0
#define RIGHT_W_VELOCITY 0.0

#define BOTTOM_NEBB_VELOCITY
#define BOTTOM_U_VELOCITY 0.0
#define BOTTOM_V_VELOCITY 0.0
#define BOTTOM_W_VELOCITY 0.0

#define TOP_NEBB_VELOCITY
#define TOP_U_VELOCITY 0.0
#define TOP_V_VELOCITY 0.0
#define TOP_W_VELOCITY 0.0

#define BACK_NEBB_VELOCITY
#define BACK_U_VELOCITY 0.0
#define BACK_V_VELOCITY 0.0
#define BACK_W_VELOCITY 0.0

#define FRONT_NEBB_VELOCITY
#define FRONT_U_VELOCITY 0.0
#define FRONT_V_VELOCITY 0.0
#define FRONT_W_VELOCITY 0.0

#define LEFT_NEBB_PRESSURE
#define LEFT_PRESSURE_RED (1.0 / 3.0)
#define LEFT_PRESSURE_BLUE 0.0

#define RIGHT_NEBB_PRESSURE
#define RIGHT_PRESSURE_RED (1.0 / 3.0)
#define RIGHT_PRESSURE_BLUE 0.0

#define BOTTOM_NEBB_PRESSURE
#define BOTTOM_PRESSURE_RED (1.0 / 3.0)
#define BOTTOM_PRESSURE_BLUE 0.0

#define TOP_NEBB_PRESSURE
#define TOP_PRESSURE_RED (1.0 / 3.0)
#define TOP_PRESSURE_BLUE 0.0

#define BACK_NEBB_PRESSURE
#define BACK_PRESSURE_RED (1.0 / 3.0)
#define BACK_PRESSURE_BLUE 0.0

#define FRONT_NEBB_PRESSURE
#define FRONT_PRESSURE_RED (1.0 / 3.0)
#define FRONT_PRESSURE_BLUE 0.0

// --------------------------
//     Initial conditions
// --------------------------
#define INI_POISEUILLE

#define INI_DROPLET
#define INI_DROPLET_R   25.0
#define INI_DROPLET_X   (0.5*(double)NX)
#define INI_DROPLET_Y   (0.5*(double)NY)
#define INI_DROPLET_Z   (0.5*(double)NZ)
#define INI_DROPLET_U   0.0
#define INI_DROPLET_V   0.0
#define INI_DROPLET_W   0.0
#define INI_DROPLET_SF  1.0

#define INI_BUBBLE
#define INI_BUBBLE_R   25.0
#define INI_BUBBLE_X   (0.5*(double)NX)
#define INI_BUBBLE_Y   (0.5*(double)NY)
#define INI_BUBBLE_Z   (0.5*(double)NZ)
#define INI_BUBBLE_U   0.0
#define INI_BUBBLE_V   0.0
#define INI_BUBBLE_W   0.0
#define INI_BUBBLE_SF  1.0

#define INI_TWOCOMPONENT_POISEUILLE
#define INI_TWOCOMPONENT_POISEUILLE_A   (0.25*(double)NY)
#define INI_TWOCOMPONENT_POISEUILLE_SF  1.0

#define INI_TWOCOMPONENT_COUETTE
#define INI_TWOCOMPONENT_COUETTE_SF     1.0

#define INI_TWODROPLETS
#define INI_TWODROPLETS_R1 40.0
#define INI_TWODROPLETS_X1 (0.25*(double)NX)
#define INI_TWODROPLETS_Y1 (0.5*(double)NY)
#define INI_TWODROPLETS_Z1 (0.5*(double)NZ)
#define INI_TWODROPLETS_U1 0.0
#define INI_TWODROPLETS_V1 0.0
#define INI_TWODROPLETS_W1 0.0
#define INI_TWODROPLETS_R2 40.0
#define INI_TWODROPLETS_X2 (0.75*(double)NX)
#define INI_TWODROPLETS_Y2 (0.5*(double)NY)
#define INI_TWODROPLETS_Z2 (0.5*(double)NZ)
#define INI_TWODROPLETS_U2 0.0
#define INI_TWODROPLETS_V2 0.0
#define INI_TWODROPLETS_W2 0.0
#define INI_TWODROPLETS_SF 1.0

#define INI_LAYERED_POISEUILLE
#define INI_LAYERED_POISEUILLE_SF  2.0

#endif
