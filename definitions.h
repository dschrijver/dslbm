#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------
//     Methods    
// ---------------
#define MRT
// #define BGK

#define COLOR_GRADIENT
// #define SHAN_CHEN

// ---------------------------
//     Boundary conditions        
// ---------------------------
#define XPERIODIC
// #define YPERIODIC
#define ZPERIODIC

// #define YPERIODIC_FLIP // For Washburn specifically

// --- Half-way Bounce-Back ---
// #define LEFT_BOUNCEBACK_NOSLIP
// #define THETA_C_LEFT 120.0
// #define XI_LEFT 0.0

// #define RIGHT_BOUNCEBACK_NOSLIP
// #define THETA_C_RIGHT 120.0
// #define XI_RIGHT 0.0

// #define BOTTOM_BOUNCEBACK_NOSLIP
// #define THETA_C_BOTTOM 120.0
// // #define XI_BOTTOM 0.0

// #define TOP_BOUNCEBACK_NOSLIP
// #define THETA_C_TOP 90.0
// // #define XI_TOP 0.0

// #define BACK_BOUNCEBACK_NOSLIP
// #define THETA_C_BACK 90.0
// #define XI_BACK 0.0

// #define FRONT_BOUNCEBACK_NOSLIP
// #define THETA_C_FRONT 90.0
// #define XI_FRONT 0.0

// #define LEFT_BOUNCEBACK_PRESSURE
// #define LEFT_PRESSURE_RED 0.0
// #define LEFT_PRESSURE_BLUE (0.4 + 3e-5)

// #define RIGHT_BOUNCEBACK_PRESSURE
// #define RIGHT_PRESSURE_RED 0.0
// #define RIGHT_PRESSURE_BLUE 0.4

// #define BOTTOM_BOUNCEBACK_PRESSURE
// #define BOTTOM_PRESSURE_RED 0.4
// #define BOTTOM_PRESSURE_BLUE 0.0

// #define TOP_BOUNCEBACK_PRESSURE
// #define TOP_PRESSURE_RED 0.0
// #define TOP_PRESSURE_BLUE 0.4

// #define BACK_BOUNCEBACK_PRESSURE
// #define BACK_PRESSURE_RED 0.4
// #define BACK_PRESSURE_BLUE 0.0

// #define FRONT_BOUNCEBACK_PRESSURE
// #define FRONT_PRESSURE_RED 0.0
// #define FRONT_PRESSURE_BLUE 0.4
// ----------------------------

// --- Non-Equilibrium Bounce-Back ---
// #define LEFT_NEBB_NOSLIP
// #define THETA_C_LEFT 120.0
// #define XI_LEFT 0.0

// #define RIGHT_NEBB_NOSLIP
// #define THETA_C_RIGHT 120.0
// #define XI_RIGHT 0.0

#define BOTTOM_NEBB_NOSLIP
#define THETA_C_BOTTOM 90.0
// #define XI_BOTTOM -0.4

#define TOP_NEBB_NOSLIP
#define THETA_C_TOP 90.0
// #define XI_TOP 0.0

// #define BACK_NEBB_NOSLIP
// #define THETA_C_BACK 90.0
// #define XI_BACK 0.0

// #define FRONT_NEBB_NOSLIP
// #define THETA_C_FRONT 90.0
// #define XI_FRONT 0.0

// #define LEFT_NEBB_PRESSURE
// #define LEFT_PRESSURE_RED 0.0
// #define LEFT_PRESSURE_BLUE (0.4 + 1e-5)

// #define RIGHT_NEBB_PRESSURE
// #define RIGHT_PRESSURE_RED 0.0
// #define RIGHT_PRESSURE_BLUE 0.4

// #define BOTTOM_NEBB_PRESSURE
// #define BOTTOM_PRESSURE_RED 0.004
// #define BOTTOM_PRESSURE_BLUE 0.0

// #define TOP_NEBB_PRESSURE
// #define TOP_PRESSURE_RED 0.0
// #define TOP_PRESSURE_BLUE 0.004

// #define BACK_NEBB_PRESSURE
// #define BACK_PRESSURE_RED 0.0
// #define BACK_PRESSURE_BLUE 0.4

// #define FRONT_NEBB_PRESSURE
// #define FRONT_PRESSURE_RED 0.0
// #define FRONT_PRESSURE_BLUE (0.4 - 1e-5)

// -----------------------------------

// --------------------------
//     Initial conditions    
// --------------------------
// #define INI_FLOATING_DROPLET
// #define R_DROPLET 15.0

// #define INI_CONTACT_ANGLE_DROPLET
// #define R_DROPLET 20.0

#define INI_POISEUILLE

// #define INI_WASHBURN
// #define H_START 25.0

// --------------------
//     Logic checks        
// --------------------
#if defined(LEFT_NEBB_NOSLIP) || defined(RIGHT_NEBB_NOSLIP) || defined(BOTTOM_NEBB_NOSLIP) || defined(TOP_NEBB_NOSLIP) || defined(BACK_NEBB_NOSLIP) || defined(FRONT_NEBB_NOSLIP) || defined(LEFT_NEBB_PRESSURE) || defined(RIGHT_NEBB_PRESSURE) || defined(BOTTOM_NEBB_PRESSURE) || defined(TOP_NEBB_PRESSURE) || defined(BACK_NEBB_PRESSURE) || defined(FRONT_NEBB_PRESSURE) || defined(RIGHT_NEBB_OUTLET)
#define WETNODE
#endif

#if defined(LEFT_BOUNCEBACK_PRESSURE) || defined(RIGHT_BOUNCEBACK_PRESSURE) || defined(BOTTOM_BOUNCEBACK_PRESSURE) || defined(TOP_BOUNCEBACK_PRESSURE) || defined(BACK_BOUNCEBACK_PRESSURE) || defined(FRONT_BOUNCEBACK_PRESSURE)
#define BOUNCEBACK_PRESSURE
#endif

#endif