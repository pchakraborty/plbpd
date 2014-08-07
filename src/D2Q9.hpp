#ifndef D2Q9_HPP
#define D2Q9_HPP

/** nVelocity **/
static int d2q9_nVelocity = 9;

/** cs2 **/
static double d2q9_cs2 = 1./3.;

/** Velocity sub-lattice */
static int d2q9_c[9][3] = { { 0, 0, 0},   //  0
			    { 1, 0, 0},   //  1
			    {-1, 0, 0},   //  2
			    { 0, 1, 0},   //  3
			    { 0,-1, 0},   //  4
			    { 1, 1, 0},   //  5
			    {-1,-1, 0},   //  6
			    { 1,-1, 0},   //  7
			    {-1, 1, 0} }; //  8

/** Weights */
static double d2q9_w[9] = { 4./9., 
			    1./9.,  1./9.,  1./9.,  1./9.,
			    1./36., 1./36., 1./36., 1./36. };

/** how many DF's to communicate */
static int d2q9_n_unk = 3;

/** which DF's to communicate in which direction */
// TODO: THIS NEEDS TO BE FIXED
static int d2q9_unk_E[3]  = {-1, -1, -1};
static int d2q9_unk_W[3]  = {-1, -1, -1};
static int d2q9_unk_N[3]  = {-1, -1, -1};
static int d2q9_unk_S[3]  = {-1, -1, -1};
static int d2q9_unk_U[3]  = {-1, -1, -1};
static int d2q9_unk_D[3]  = {-1, -1, -1};

/** reverse - for noslip boundary condition */
static int d2q9_reverse[19] = { 0,
				2,  1,
				4,  3,
				6,  5,
				8,  7 };

#endif
