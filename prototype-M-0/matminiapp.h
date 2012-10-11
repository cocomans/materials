/*CoCoMANS miniapp
version 1.0*/

//common functions and variables
inline void* xmalloc(size_t size, const char *what) {
  void* out;
  out = malloc(size);
  if(size && !out) {
    fprintf(stderr,"null pointer in malloc for %s\n",what);
    abort();
  }
  return out;
}

inline void* xrealloc(void *ptr, size_t size, const char *what) {
  void* out;
  out = realloc(ptr, size);
  if(size && !out) {
    fprintf(stderr,"null pointer in realloc for %s\n",what);
    abort();
  }
  return out;
}

/*index function for pointers used instead of arrays
array[lmax][kmax][jmax][imax]
is declared as *array of size imax*jmax*kmax*lmax
accessing array[l][k][j][i] is equivalent to using
array[index(i,imax,j,jmax,k,kmax,l)] or *(array+index(i,imax,j,jmax,k,kmax,l))*/
inline int index(int i,int imax,int j,int jmax = 0,int k = 0,int kmax = 0, int l = 0) {
  return i + imax*(j + jmax*(k + kmax*l));
}

//-----------------------
//Common variables
//-----------------------

double *sigma_11;

//-----------------------
//MPM variables
//-----------------------
int ncells, parsPerCell, nnodes, np_mpm;
double loadingTime = 0.050;
double length_mpm, dt_mpm, dx_mpm, parSize;
double YoungM;
double *nodex, *nodeVel, *nodeMass;
double *xp, *mp_mpm, *velp, *langVel, *grdVol, *lp,*l0_mpm,*strain;
double parVol, gravity, time_mpm;
double targetG = 0.0;
double timeOut; 
double sinceLastOut = 0;
double iclast = -0.1; //initial condition on last node in m/s
double icfirst = 0; //initial condition on first node in m/s
double bclast = -0.1; //boundary condition on last node in m/s
double bcfirst = 0; //boundary condition on first node in m/s

//-----------------------
//MD variables
//-----------------------

/* Initialized variables*/
const int mp_md = 2504;     // Should be a multiple of 8 but *not* of 4096 (32KB cache/8 ways per cache set)

template <size_t NUM_ELTS>
class XYZvector {
 public:
  double xyz[3][NUM_ELTS];
};

/*Simluation variabes*/
double time_md,dt_md,maxtime;
double xyzL[3],xyzL0[3];

int np_md,it,nstep,relaxation,par;

int freq;

int nsecout;
double vx_Mean;
double penergy,kenergy,strain_md,final_strain;

/*Atom related arrays*/
XYZvector<mp_md> pos;
XYZvector<mp_md> previous;

XYZvector<mp_md> velocity;
XYZvector<mp_md> velocity_h;   // Velocity at half timesteps

double rmass;

XYZvector<mp_md> force;

double gamxx0,epsilon,sigma;

/*Arrays needed for stress calculations*/
double vgrad[3][3];
double trelax;
double stress[3][3];

// We promise not to use fewer than this many <whatever>s.
#define MIN_NP_MD 1024

