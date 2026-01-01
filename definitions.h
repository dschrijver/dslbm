#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions        
// ---------------------------
// #define XPERIODIC
// #define YPERIODIC
// #define ZPERIODIC

#define LEFT_BOUNCEBACK
#define THETA_C_LEFT 90.0

#define RIGHT_BOUNCEBACK
#define THETA_C_RIGHT 90.0

// #define BOTTOM_BOUNCEBACK
// #define THETA_C_BOTTOM 90.0

// #define TOP_BOUNCEBACK
// #define THETA_C_TOP 90.0

#define BACK_BOUNCEBACK
#define THETA_C_BACK 90.0

#define FRONT_BOUNCEBACK
#define THETA_C_FRONT 90.0

#define BOTTOM_NEBB_NOSLIP
#define THETA_C_BOTTOM 90.0

#define TOP_NEBB_NOSLIP
#define THETA_C_TOP 90.0

// --------------------------
//     Initial conditions    
// --------------------------
// #define INI_FLOATING_DROPLET
// #define R_DROPLET 15.0

#define INI_CONTACT_ANGLE_DROPLET
#define R_DROPLET 15.0

// #define INI_POISEUILLE

#endif