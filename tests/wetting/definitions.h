#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define XPERIODIC
#define ZPERIODIC

// --- Wetting ---
#define THETA_C 150.0

#define BOTTOM_HWBB_NOSLIP
#define TOP_HWBB_NOSLIP

// --------------------------
//     Initial conditions
// --------------------------
#define INI_DROPLET
#define INI_DROPLET_R   40.0
#define INI_DROPLET_X   (0.5*(double)NX)
#define INI_DROPLET_Y   40.0
#define INI_DROPLET_Z   (0.5*(double)NZ)
#define INI_DROPLET_U   0.0
#define INI_DROPLET_V   0.0
#define INI_DROPLET_W   0.0
#define INI_DROPLET_SF  2.0

#endif
