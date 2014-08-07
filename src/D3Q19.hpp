#ifndef D3Q19_HPP
#define D3Q19_HPP

/** nVelocity **/
static int d3q19_nVelocity = 19;

/** cs2 **/
static double d3q19_cs2 = 1./3.;

/** Velocity sub-lattice */
static int d3q19_c[19][3] = { { 0, 0, 0},   //  0
		              { 1, 0, 0},   //  1
			      {-1, 0, 0},   //  2
			      { 0, 1, 0},   //  3
			      { 0,-1, 0},   //  4
			      { 0, 0, 1},   //  5
			      { 0, 0,-1},   //  6
			      { 1, 1, 0},   //  7
			      {-1,-1, 0},   //  8
			      { 1,-1, 0},   //  9
			      {-1, 1, 0},   // 10
			      { 1, 0, 1},   // 11
			      {-1, 0,-1},   // 12
			      { 1, 0,-1},   // 13
			      {-1, 0, 1},   // 14
			      { 0, 1, 1},   // 15
			      { 0,-1,-1},   // 16
			      { 0, 1,-1},   // 17
			      { 0,-1, 1} }; // 18

/** Weights */
static double d3q19_w[19] = { 1./3., 
			      1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
			      1./36., 1./36., 1./36., 1./36., 
			      1./36., 1./36., 1./36., 1./36.,
			      1./36., 1./36., 1./36., 1./36. };

/** how many DF's to communicate */
static int d3q19_n_unk = 5;

/** which DF's to communicate in which direction */
static int d3q19_unk_E[5]  = {2,  8, 10, 12, 14};
static int d3q19_unk_W[5]  = {1,  7,  9, 11, 13};
static int d3q19_unk_N[5]  = {4,  8,  9, 16, 18};
static int d3q19_unk_S[5]  = {3,  7, 10, 15, 17};
static int d3q19_unk_U[5]  = {6, 12, 13, 16, 17};
static int d3q19_unk_D[5]  = {5, 11, 14, 15, 18};

/** reverse - for noslip boundary condition */
static int d3q19_reverse[19] = { 0,
				 2,  1,
				 4,  3,
				 6,  5,
				 8,  7,
				10,  9,
				12, 11,
				14, 13,
				16, 15,
				18, 17};

#endif
