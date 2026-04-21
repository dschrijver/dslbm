#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define XPERIODIC
#define ZPERIODIC

// --- Wetting ---
#define THETA_C 160.0

#define BOTTOM_REGNEE_VELOCITY
#define BOTTOM_U_VELOCITY 0.0
#define BOTTOM_V_VELOCITY 0.0
#define BOTTOM_W_VELOCITY 0.0

#define TOP_REGNEE_VELOCITY
#define TOP_U_VELOCITY 0.0
#define TOP_V_VELOCITY 0.0
#define TOP_W_VELOCITY 0.0

// --------------------------
//     Initial conditions
// --------------------------
#define INI_DROPLET
#define INI_DROPLET_R   25.0
#define INI_DROPLET_X   (0.5*(double)NX)
#define INI_DROPLET_Y   0.0
#define INI_DROPLET_Z   (0.5*(double)NZ)
#define INI_DROPLET_U   0.0
#define INI_DROPLET_V   0.0
#define INI_DROPLET_W   0.0
#define INI_DROPLET_SF  2.0

#endif
