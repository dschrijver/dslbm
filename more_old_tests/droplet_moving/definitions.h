#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions
// ---------------------------
#define XPERIODIC
#define YPERIODIC
#define ZPERIODIC

// --------------------------
//     Initial conditions
// --------------------------
#define INI_DROPLET
#define INI_DROPLET_R   15.0
#define INI_DROPLET_X   (0.5*(double)NX)
#define INI_DROPLET_Y   (0.5*(double)NY)
#define INI_DROPLET_Z   (0.5*(double)NZ)
#define INI_DROPLET_U   0.0
#define INI_DROPLET_V   0.0
#define INI_DROPLET_W   0.0
#define INI_DROPLET_SF  2.0

#endif
