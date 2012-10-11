#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "MDminiapp.h"
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <algorithm>

//-----------------
//MD
//-----------------

//define a method to shift arg by amount of n*prd
//and into (-prd/2,prd/2].
double lShift(double arg, double prd) {
  double arg_after;
  arg_after = arg + floor(0.5 - arg / prd) * prd;
  return arg_after;
}

// Set some variables to default value
void setDefault() {

  // initialize the macroscopic velocity gradient.
#ifdef __INTEL_COMPILER
  vgrad[:][:] = 0.0;
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      vgrad[i][j] = 0.0;
#endif
}

   // Get the input parameters from the input file
void dataInput() {

  char dummy[50];
  // open a input file
  FILE *in = fopen("dry_flow.dat","r");
  // start to read data
  fscanf(in,"%s",dummy);
  fscanf(in,"%lf\n",&epsilon);
  //fprintf(stdout,"%s = %f\n",dummy,epsilon);
  fscanf(in,"%s",dummy);
  fscanf(in,"%lf",&sigma);
  //fprintf(stdout,"%s = %f\n",dummy,sigma);
  fscanf(in,"%s",dummy);
  fscanf(in,"%lf",&dt_md);
  //fprintf(stdout,"%s = %f\n",dummy,dt_md);
  fscanf(in,"%s",dummy);
  fscanf(in,"%lf",&trelax);
  //fprintf(stdout,"%s = %f\n",dummy,trelax);
  fscanf(in,"%s",dummy);
  fscanf(in,"%lf",&vx_Mean);
  //fprintf(stdout,"%s = %f\n",dummy,vx_Mean);
  fscanf(in,"%s",dummy);
  fscanf(in,"%d",&freq);
  //fprintf(stdout,"%s = %d\n",dummy,freq);
  fscanf(in,"%s",dummy);
  fscanf(in,"%lf",&maxtime);
  //fprintf(stdout,"%s = %f\n",dummy,maxtime);

  // close the file
  fclose(in);
}

void setInit() {

  int atomid;
  FILE *in;
  // read the initial conditions from the initial file.
  if (relaxation == 1)
    in = fopen("pvm.dat","r");
  else 
    in = fopen("restart.dat","r");
  // read the first line of the file
  fscanf(in,"%d",&np_md);

  if (np_md > mp_md) {
    printf("maximum number of particles allowed exceeded\n");
    abort();
  }

  fscanf(in,"%lf %lf %lf",&xyzL[0],&xyzL[1],&xyzL[2]);

  xyzL0[0] = xyzL[0];
  xyzL0[1] = xyzL[1];
  xyzL0[2] = xyzL[2];

  strain_md = 0;
  fscanf(in,"%lf",&rmass);
   // then read the rest of the data.
  // the inital positions and radius of each particle.
  for (int i = 0; i < np_md; i++)
    fscanf(in,"%d %lf %lf %lf",&atomid,&pos.xyz[0][i],&pos.xyz[1][i],&pos.xyz[2][i]);

  // the initial velocities.
  for (int i = 0; i < np_md; i++)
    fscanf(in,"%d %lf %lf %lf",&atomid,&velocity.xyz[0][i],&velocity.xyz[1][i],&velocity.xyz[2][i]);

  // close the file
  fclose(in);

  nsecout = 0;
}

   // this method update the macroscopic velocity gradient
void newVgrad() {

   if (relaxation == 1) {
    vgrad[0][0] = 0.0;
  }
  else {
    if (strain_md < final_strain) {
      // initial strain rate
      vgrad[0][0] = gamxx0;
    }
    else {
      vgrad[0][0] = 0.0;
    }
  }
}

// this method will adjust the particle velocities according to the
// macroscopic velocity gradient using the least squares method to fit
// to a straight line.
void velAdjust() {

  double grd[3][3];
  double vmean[3];
  double sxx, sxy, sxz, syy, syz, szz;
  double xvx, yvx, zvx;
  double xvy, yvy, zvy;
  double xvz, yvz, zvz;
  double g[3][3];
  double det;

#ifdef __INTEL_COMPILER
  // calculate the mean velocity in each direction
# pragma loop count min(MIN_NP_MD)
  vmean[0] = __sec_reduce_add(velocity.xyz[0][0:np_md]) / np_md;
  vmean[1] = __sec_reduce_add(velocity.xyz[1][0:np_md]) / np_md;
  vmean[2] = __sec_reduce_add(velocity.xyz[2][0:np_md]) / np_md;

  // subtract the mean
# pragma loop count min(MIN_NP_MD)
  velocity.xyz[0][0:np_md] -= vmean[0];
  velocity.xyz[1][0:np_md] -= vmean[1];
  velocity.xyz[2][0:np_md] -= vmean[2];

# pragma loop count min(MIN_NP_MD)
  sxx = __sec_reduce_add(pos.xyz[0][0:np_md] * pos.xyz[0][0:np_md]);
  sxy = __sec_reduce_add(pos.xyz[0][0:np_md] * pos.xyz[1][0:np_md]);
  sxz = __sec_reduce_add(pos.xyz[0][0:np_md] * pos.xyz[2][0:np_md]);
  syy = __sec_reduce_add(pos.xyz[1][0:np_md] * pos.xyz[1][0:np_md]);
  syz = __sec_reduce_add(pos.xyz[1][0:np_md] * pos.xyz[2][0:np_md]);
  szz = __sec_reduce_add(pos.xyz[2][0:np_md] * pos.xyz[2][0:np_md]);

  xvx = __sec_reduce_add(pos.xyz[0][0:np_md] * velocity.xyz[0][0:np_md]);
  yvx = __sec_reduce_add(pos.xyz[1][0:np_md] * velocity.xyz[0][0:np_md]);
  zvx = __sec_reduce_add(pos.xyz[2][0:np_md] * velocity.xyz[0][0:np_md]);

  xvy = __sec_reduce_add(pos.xyz[0][0:np_md] * velocity.xyz[1][0:np_md]);
  yvy = __sec_reduce_add(pos.xyz[1][0:np_md] * velocity.xyz[1][0:np_md]);
  zvy = __sec_reduce_add(pos.xyz[2][0:np_md] * velocity.xyz[1][0:np_md]);

  xvz = __sec_reduce_add(pos.xyz[0][0:np_md] * velocity.xyz[2][0:np_md]);
  yvz = __sec_reduce_add(pos.xyz[1][0:np_md] * velocity.xyz[2][0:np_md]);
  zvz = __sec_reduce_add(pos.xyz[2][0:np_md] * velocity.xyz[2][0:np_md]);
#else
  vmean[0] = 0.0;
  vmean[1] = 0.0;
  vmean[2] = 0.0;

  // calculate the mean velocity in each direction
  double vxmean = 0.0;
  double vymean = 0.0;
  double vzmean = 0.0;
# pragma omp parallel for schedule(static) reduction(+:vxmean,vymean,vzmean)
# pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    vxmean += velocity.xyz[0][i];
    vymean += velocity.xyz[1][i];            //vector vmean + vector velocities
    vzmean += velocity.xyz[2][i];
  }
  vmean[0] = vxmean / np_md;
  vmean[1] = vymean / np_md;             //vector vmean / scalar np_md
  vmean[2] = vzmean / np_md;

  // subtract the mean
# pragma omp parallel for schedule(static)
# pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    velocity.xyz[0][i] -= vmean[0];
    velocity.xyz[1][i] -= vmean[1];           //vector velocities - vector vmean
    velocity.xyz[2][i] -= vmean[2];
  }

  sxx = 0.0;
  sxy = 0.0;
  sxz = 0.0;
  syy = 0.0;
  syz = 0.0;
  szz = 0.0;

  xvx = 0.0;
  yvx = 0.0;
  zvx = 0.0;

  xvy = 0.0;
  yvy = 0.0;
  zvy = 0.0;

  xvz = 0.0;
  yvz = 0.0;
  zvz = 0.0;

# pragma omp parallel for schedule(static) reduction(+:sxx,sxy,sxz,syy,syz,szz,xvx,yvx,zvx,xvy,yvy,zvy,xvz,yvz,zvz)
# pragma loop count min(MIN_NP_MD)
  for (int ip = 0; ip < np_md; ip++) {
    sxx += pos.xyz[0][ip] * pos.xyz[0][ip];
    sxy += pos.xyz[0][ip] * pos.xyz[1][ip];
    sxz += pos.xyz[0][ip] * pos.xyz[2][ip];          //matrix S = vector (column) position * vector (row) position (=transpose of vector (column position))
    syy += pos.xyz[1][ip] * pos.xyz[1][ip];
    syz += pos.xyz[1][ip] * pos.xyz[2][ip];
    szz += pos.xyz[2][ip] * pos.xyz[2][ip];

    xvx += pos.xyz[0][ip] * velocity.xyz[0][ip];
    yvx += pos.xyz[1][ip] * velocity.xyz[0][ip];        //matrix xvx += vector (column) velocities * vector (row) positions
    zvx += pos.xyz[2][ip] * velocity.xyz[0][ip];

    xvy += pos.xyz[0][ip] * velocity.xyz[1][ip];
    yvy += pos.xyz[1][ip] * velocity.xyz[1][ip];
    zvy += pos.xyz[2][ip] * velocity.xyz[1][ip];

    xvz += pos.xyz[0][ip] * velocity.xyz[2][ip];
    yvz += pos.xyz[1][ip] * velocity.xyz[2][ip];
    zvz += pos.xyz[2][ip] * velocity.xyz[2][ip];
  }
#endif

  det = sxx * syy * szz + sxy * syz * sxz + sxz * syz * sxy -          //determinant of S
    sxz * syy * sxz - syz * syz * sxx - szz * sxy * sxy;

  g[0][0] = xvx * syy * szz + yvx * syz * sxz + zvx * syz * sxy -
    sxz * syy * zvx - syz * syz * xvx - szz * yvx * sxy;

  g[0][1] = sxx * yvx * szz + sxy * zvx * sxz + sxz * syz * xvx -
    sxz * yvx * sxz - syz * zvx * sxx - szz * sxy * xvx;

  g[0][2] = sxx * syy * zvx + sxy * syz * xvx + sxz * yvx * sxy -    //g[i][j] is the determinant of the matrix S where row j has been replaced
    xvx * syy * sxz - yvx * syz * sxx - zvx * sxy * sxy;             //by row i of matrix xvx

  g[1][0] = xvy * syy * szz + yvy * syz * sxz + zvy * syz * sxy -
    sxz * syy * zvy - syz * syz * xvy - szz * yvy * sxy;

  g[1][1] = sxx * yvy * szz + sxy * zvy * sxz + sxz * syz * xvy -
    sxz * yvy * sxz - syz * zvy * sxx - szz * sxy * xvy;

  g[1][2] = sxx * syy * zvy + sxy * syz * xvy + sxz * yvy * sxy -
    xvy * syy * sxz - yvy * syz * sxx - zvy * sxy * sxy;

  g[2][0] = xvz * syy * szz + yvz * syz * sxz + zvz * syz * sxy -
    sxz * syy * zvz - syz * syz * xvz - szz * yvz * sxy;

  g[2][1] = sxx * yvz * szz + sxy * zvz * sxz + sxz * syz * xvz -
    sxz * yvz * sxz - syz * zvz * sxx - szz * sxy * xvz;

  g[2][2] = sxx * syy * zvz + sxy * syz * xvz + sxz * yvz * sxy -
    xvz * syy * sxz - yvz * syz * sxx - zvz * sxy * sxy;

#ifdef __INTEL_COMPILER
  g[:][:] /= det;
#else
  g[0][0] /= det;
  g[0][1] /= det;
  g[0][2] /= det;

  g[1][0] /= det;          //matrix g / scalar det
  g[1][1] /= det;
  g[1][2] /= det;

  g[2][0] /= det;
  g[2][1] /= det;
  g[2][2] /= det;
#endif

  // update the new velocity gradient
#ifdef __INTEL_COMPILER
  grd[:][:] = vgrad[:][:] - g[:][:];
#else
  grd[0][0] = vgrad[0][0] - g[0][0];
  grd[0][1] = vgrad[0][1] - g[0][1];
  grd[0][2] = vgrad[0][2] - g[0][2];

  grd[1][0] = vgrad[1][0] - g[1][0];
  grd[1][1] = vgrad[1][1] - g[1][1];       //matrix grd = matrix vgrad - matrix g
  grd[1][2] = vgrad[1][2] - g[1][2];

  grd[2][0] = vgrad[2][0] - g[2][0];
  grd[2][1] = vgrad[2][1] - g[2][1];
  grd[2][2] = vgrad[2][2] - g[2][2];
#endif

#ifdef __INTEL_COMPILER
  // adjust the velocities
  for (int dim = 0; dim < 3; dim++)
# pragma loop count min(MIN_NP_MD)
    velocity.xyz[dim][0:np_md] += grd[dim][0] * pos.xyz[0][0:np_md] + grd[dim][1] * pos.xyz[1][0:np_md] + grd[dim][2] * pos.xyz[2][0:np_md];

  // adjust the velocities according the mean velocities
  velocity.xyz[0][0:np_md] += vx_Mean;

  // now adjust the half-step velocity
  // and calculate the mean velocity in each direction
  vmean[0] = __sec_reduce_add(velocity_h.xyz[0][0:np_md]) / np_md;
  vmean[1] = __sec_reduce_add(velocity_h.xyz[1][0:np_md]) / np_md;
  vmean[2] = __sec_reduce_add(velocity_h.xyz[2][0:np_md]) / np_md;

  // subtract the mean
# pragma loop count min(MIN_NP_MD)
  velocity_h.xyz[0][0:np_md] -= vmean[0];
  velocity_h.xyz[1][0:np_md] -= vmean[1];
  velocity_h.xyz[2][0:np_md] -= vmean[2];

# pragma loop count min(MIN_NP_MD)
  sxx = __sec_reduce_add(pos.xyz[0][0:np_md] * pos.xyz[0][0:np_md]);
  sxy = __sec_reduce_add(pos.xyz[0][0:np_md] * pos.xyz[1][0:np_md]);
  sxz = __sec_reduce_add(pos.xyz[0][0:np_md] * pos.xyz[2][0:np_md]);
  syy = __sec_reduce_add(pos.xyz[1][0:np_md] * pos.xyz[1][0:np_md]);
  syz = __sec_reduce_add(pos.xyz[1][0:np_md] * pos.xyz[2][0:np_md]);
  szz = __sec_reduce_add(pos.xyz[2][0:np_md] * pos.xyz[2][0:np_md]);

  xvx = __sec_reduce_add(pos.xyz[0][0:np_md] * velocity_h.xyz[0][0:np_md]);
  yvx = __sec_reduce_add(pos.xyz[1][0:np_md] * velocity_h.xyz[0][0:np_md]);
  zvx = __sec_reduce_add(pos.xyz[2][0:np_md] * velocity_h.xyz[0][0:np_md]);

  xvy = __sec_reduce_add(pos.xyz[0][0:np_md] * velocity_h.xyz[1][0:np_md]);
  yvy = __sec_reduce_add(pos.xyz[1][0:np_md] * velocity_h.xyz[1][0:np_md]);
  zvy = __sec_reduce_add(pos.xyz[2][0:np_md] * velocity_h.xyz[1][0:np_md]);

  xvz = __sec_reduce_add(pos.xyz[0][0:np_md] * velocity_h.xyz[2][0:np_md]);
  yvz = __sec_reduce_add(pos.xyz[1][0:np_md] * velocity_h.xyz[2][0:np_md]);
  zvz = __sec_reduce_add(pos.xyz[2][0:np_md] * velocity_h.xyz[2][0:np_md]);
#else
# pragma omp parallel for schedule(static)
# pragma loop count min(MIN_NP_MD)
  // adjust the velocities
  for (int ip = 0; ip < np_md; ip++) {
    velocity.xyz[0][ip] += grd[0][0] * pos.xyz[0][ip] + grd[0][1] * pos.xyz[1][ip] + grd[0][2] * pos.xyz[2][ip]; //vector velocities += matrix grd * vector (column) positions
    velocity.xyz[1][ip] += grd[1][0] * pos.xyz[0][ip] + grd[1][1] * pos.xyz[1][ip] + grd[1][2] * pos.xyz[2][ip];
    velocity.xyz[2][ip] += grd[2][0] * pos.xyz[0][ip] + grd[2][1] * pos.xyz[1][ip] + grd[2][2] * pos.xyz[2][ip];
  }

  // adjust the velocities according the mean velocities
# pragma omp parallel for schedule(static)
# pragma loop count min(MIN_NP_MD)
  for (int ip = 0; ip < np_md; ip++)
    velocity.xyz[0][ip] += vx_Mean;

  // now adjust the half-step velocity
  vmean[0] = 0.0;
  vmean[1] = 0.0;
  vmean[2] = 0.0;

  // calculate the mean velocity in each direction
  vxmean = 0.0;
  vymean = 0.0;
  vzmean = 0.0;
# pragma omp parallel for schedule(static) reduction(+:vxmean,vymean,vzmean)
# pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    vxmean += velocity.xyz[0][i];
    vymean += velocity.xyz[1][i];            //vector vmean + vector velocities
    vzmean += velocity.xyz[2][i];
  }
  vmean[0] = vxmean / np_md;
  vmean[1] = vymean / np_md;             //vector vmean / scalar np_md
  vmean[2] = vzmean / np_md;

  // subtract the mean
# pragma omp parallel for schedule(static)
# pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    velocity_h.xyz[0][i] -= vmean[0];
    velocity_h.xyz[1][i] -= vmean[1];
    velocity_h.xyz[2][i] -= vmean[2];
  }

  sxx = 0.0;
  sxy = 0.0;
  sxz = 0.0;
  syy = 0.0;
  syz = 0.0;
  szz = 0.0;

  xvx = 0.0;
  yvx = 0.0;
  zvx = 0.0;

  xvy = 0.0;
  yvy = 0.0;
  zvy = 0.0;

  xvz = 0.0;
  yvz = 0.0;
  zvz = 0.0;

# pragma omp parallel for schedule(static) reduction(+:sxx,sxy,sxz,syy,syz,szz,xvx,yvx,zvx,xvy,yvy,zvy,xvz,yvz,zvz)
# pragma loop count min(MIN_NP_MD)
  for (int ip = 0; ip < np_md; ip++) {
    sxx += pos.xyz[0][ip] * pos.xyz[0][ip];
    sxy += pos.xyz[0][ip] * pos.xyz[1][ip];
    sxz += pos.xyz[0][ip] * pos.xyz[2][ip];
    syy += pos.xyz[1][ip] * pos.xyz[1][ip];
    syz += pos.xyz[1][ip] * pos.xyz[2][ip];
    szz += pos.xyz[2][ip] * pos.xyz[2][ip];

    xvx += pos.xyz[0][ip] * velocity_h.xyz[0][ip];
    yvx += pos.xyz[1][ip] * velocity_h.xyz[0][ip];
    zvx += pos.xyz[2][ip] * velocity_h.xyz[0][ip];

    xvy += pos.xyz[0][ip] * velocity_h.xyz[1][ip];
    yvy += pos.xyz[1][ip] * velocity_h.xyz[1][ip];
    zvy += pos.xyz[2][ip] * velocity_h.xyz[1][ip];

    xvz += pos.xyz[0][ip] * velocity_h.xyz[2][ip];
    yvz += pos.xyz[1][ip] * velocity_h.xyz[2][ip];
    zvz += pos.xyz[2][ip] * velocity_h.xyz[2][ip];
  }
#endif

  det = sxx * syy * szz + sxy * syz * sxz + sxz * syz * sxy -
    sxz * syy * sxz - syz * syz * sxx - szz * sxy * sxy;

  g[0][0] = xvx * syy * szz + yvx * syz * sxz + zvx * syz * sxy -
    sxz * syy * zvx - syz * syz * xvx - szz * yvx * sxy;

  g[0][1] = sxx * yvx * szz + sxy * zvx * sxz + sxz * syz * xvx -
    sxz * yvx * sxz - syz * zvx * sxx - szz * sxy * xvx;

  g[0][2] = sxx * syy * zvx + sxy * syz * xvx + sxz * yvx * sxy -
    xvx * syy * sxz - yvx * syz * sxx - zvx * sxy * sxy;

  g[1][0] = xvy * syy * szz + yvy * syz * sxz + zvy * syz * sxy -
    sxz * syy * zvy - syz * syz * xvy - szz * yvy * sxy;

  g[1][1] = sxx * yvy * szz + sxy * zvy * sxz + sxz * syz * xvy -
    sxz * yvy * sxz - syz * zvy * sxx - szz * sxy * xvy;

  g[1][2] = sxx * syy * zvy + sxy * syz * xvy + sxz * yvy * sxy -
    xvy * syy * sxz - yvy * syz * sxx - zvy * sxy * sxy;

  g[2][0] = xvz * syy * szz + yvz * syz * sxz + zvz * syz * sxy -
    sxz * syy * zvz - syz * syz * xvz - szz * yvz * sxy;

  g[2][1] = sxx * yvz * szz + sxy * zvz * sxz + sxz * syz * xvz -
    sxz * yvz * sxz - syz * zvz * sxx - szz * sxy * xvz;

  g[2][2] = sxx * syy * zvz + sxy * syz * xvz + sxz * yvz * sxy -
    xvz * syy * sxz - yvz * syz * sxx - zvz * sxy * sxy;

#ifdef __INTEL_COMPILER
  // update the new velocity gradient
  g[:][:] /= det;
  grd[:][:] = vgrad[:][:] - g[:][:];

  // adjust the velocities according the mean velocities
  for (int dim = 0; dim < 3; dim++)
    velocity_h.xyz[dim][0:np_md] += grd[dim][0] * pos.xyz[0][0:np_md] + grd[dim][1] * pos.xyz[1][0:np_md] + grd[dim][2] * pos.xyz[2][0:np_md];
  velocity_h.xyz[0][0:np_md] += vx_Mean;
#else
  // update the new velocity gradient
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 3; c++) {
      g[r][c] /= det;
      grd[r][c] = vgrad[r][c] - g[r][c];
    }

  // adjust the velocities according the mean velocities
# pragma omp parallel for schedule(static)
# pragma loop count min(MIN_NP_MD)
  for (int ip = 0; ip < np_md; ip++) {
    for (int dim = 0; dim < 3; dim++)
      velocity_h.xyz[dim][ip] += grd[dim][0] * pos.xyz[0][ip] + grd[dim][1] * pos.xyz[1][ip] + grd[dim][2] * pos.xyz[2][ip];
    velocity_h.xyz[0][ip] += vx_Mean;
  }
#endif
}

// define a method to calculate the particle forces
void parForces() {

#ifdef __INTEL_COMPILER
# pragma loop count min(MIN_NP_MD)
  force.xyz[:][0:np_md] = 0.0;
#else
# pragma omp parallel for schedule(static)
# pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    force.xyz[0][i] = 0.0;
    force.xyz[1][i] = 0.0;
    force.xyz[2][i] = 0.0;
  }
#endif

  double xyzL_local[3] = {xyzL[0], xyzL[1], xyzL[2]};

#pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    double force_prod_i[3][np_md];
    for (int j = 0; j < i; j++) {
      // normal force
      double dis[3] = {
        (pos.xyz[0][i] - pos.xyz[0][j]),
        (pos.xyz[1][i] - pos.xyz[1][j]),
        (pos.xyz[2][i] - pos.xyz[2][j])
      };
      if (fabs(dis[0]) > 0.5*xyzL_local[0]) {
        double s = 1.0;
        if (dis[0] < 0)
          s = -1.0;
        dis[0] -= s*xyzL_local[0];
      }
      if (fabs(dis[1]) > 0.5*xyzL_local[1]) {
        double s = 1.0;
        if (dis[1] < 0)
          s = -1.0;
        dis[1] -= s*xyzL_local[1];
      }
      if (fabs(dis[2]) > 0.5*xyzL_local[2]) {
        double s = 1.0;
        if (dis[2] < 0)
          s = -1.0;
        dis[2] -= s*xyzL_local[2];
      }

      float distance2 = dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2];
      float inv_distance2 = 1.0f/distance2;
      double sigma6 = pow(sigma, 6.0);
      double sigma12 = sigma6*sigma6;
      double inv_distance8 = pow(inv_distance2, 4);
      double inv_distance14 = pow(inv_distance2, 7);
      double fnorm = 24. * epsilon * (2.0*sigma12 * inv_distance14 - sigma6 * inv_distance8);

      force_prod_i[0][j] = fnorm * dis[0];
      force_prod_i[1][j] = fnorm * dis[1];         //vector force += scalar fnorm * vector dis
      force_prod_i[2][j] = fnorm * dis[2];
      force.xyz[0][j] -= fnorm * dis[0];
      force.xyz[1][j] -= fnorm * dis[1];         //vector force += scalar fnorm * vector dis
      force.xyz[2][j] -= fnorm * dis[2];
    }

    for (int j = 0; j < i; j++) {
      force.xyz[0][i] += force_prod_i[0][j];
      force.xyz[1][i] += force_prod_i[1][j];
      force.xyz[2][i] += force_prod_i[2][j];
    }
  }
}

void energy() {
  penergy = 0.;
  kenergy = 0.;

  //energy calculation
#pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    //fprintf(stdout,"atom %d coordinates: %f %f %f\n",i,coord[i].x,coord[i].y,coord[i].z);
    //fprintf(stdout,"atom %d velocity: %f %f %f\n",i,v[i].x,v[i].y,v[i].z);
    for (int j = 0; j < i; j++) {
      double dis[3] = {
        (pos.xyz[0][j] - pos.xyz[0][i]),
        (pos.xyz[1][j] - pos.xyz[1][i]),
        (pos.xyz[2][j] - pos.xyz[2][i])
      };
      if(abs(dis[0]) > 0.5*xyzL[0]) dis[0] = -copysign(xyzL[0],dis[0]) + dis[0];   //periodicity
      if(abs(dis[1]) > 0.5*xyzL[1]) dis[1] = -copysign(xyzL[1],dis[1]) + dis[1];
      if(abs(dis[2]) > 0.5*xyzL[2]) dis[2] = -copysign(xyzL[2],dis[2]) + dis[2];
      double distance = sqrt(dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]);
      
      double pow6 = pow(sigma/distance, 6);
      penergy += 4. * epsilon * (pow6*pow6 - pow6);
    }
#ifndef __INTEL_COMPILER
    kenergy += 0.5 * rmass * (velocity.xyz[0][i]*velocity.xyz[0][i] + velocity.xyz[1][i]*velocity.xyz[1][i] + velocity.xyz[2][i]*velocity.xyz[2][i]);          //kenergy += 0.5 * scalar rmass * dotproduct(velocities,velocities)
#endif
  }
#ifdef __INTEL_COMPILER
# pragma loop count min(MIN_NP_MD)
  kenergy = 0.5 * rmass * __sec_reduce_add(velocity.xyz[0][0:np_md]*velocity.xyz[0][0:np_md] + velocity.xyz[1][0:np_md]*velocity.xyz[1][0:np_md] +velocity.xyz[2][0:np_md]*velocity.xyz[2][0:np_md]);
#endif
}

// this method is the integration of particles using the Leapfrog method.
void parMove() {

  //at the very begin of simulation, use Euler algorithm to get first half-step velocity,
  //and Leapfrog for the first positions.
  if (time_md <= 0.0) {
    parForces();
#ifdef __INTEL_COMPILER
# pragma loop count min(MIN_NP_MD)
    velocity_h.xyz[:][0:np_md] = velocity.xyz[:][0:np_md] + force.xyz[:][0:np_md] * dt_md / (2. * rmass);
    previous.xyz[:][0:np_md] = pos.xyz[:][0:np_md];
    pos.xyz[:][0:np_md] += velocity_h.xyz[:][0:np_md] * dt_md;
#else
# pragma omp parallel for collapse(2) schedule(static)
# pragma loop count min(MIN_NP_MD)
    for (int i = 0; i < np_md; i++)
      for (int dim = 0; dim < 3; dim++) {
        velocity_h.xyz[dim][i] = velocity.xyz[dim][i] + force.xyz[dim][i] * dt_md / (2. * rmass);      //vector vh = vector velocities + vector forces * scalar (dt_md / (2*rmass))
        previous.xyz[dim][i] = pos.xyz[dim][i];             //vector previous = vector positions
        pos.xyz[dim][i] += velocity_h.xyz[dim][i] * dt_md;      // vector positions += vector vh * scalar dt_md
      }
#endif
  }
  else {
    // other than first step, use Leapfrog method to determine the force.
    //!!!! velocities at integer timesteps lag behind the positions by dt_md
    parForces();
#ifdef __INTEL_COMPILER
# pragma loop count min(MIN_NP_MD)
    velocity.xyz[:][0:np_md] = velocity_h.xyz[:][0:np_md] / 2.;
    velocity_h.xyz[:][0:np_md] += force.xyz[:][0:np_md] * dt_md / rmass;   //at time t+0.5dt_md
    velocity.xyz[:][0:np_md] += velocity_h.xyz[:][0:np_md] / 2.;      //at time t
    previous.xyz[:][0:np_md] = pos.xyz[:][0:np_md];               //at time t
    pos.xyz[:][0:np_md] += velocity_h.xyz[:][0:np_md] * dt_md;            //at time t+dt_md
#else
# pragma omp parallel for collapse(2) schedule(static)
# pragma loop count min(MIN_NP_MD)
    for (int i = 0; i < np_md; i++)
      for (int dim = 0; dim < 3; dim++) {
        velocity.xyz[dim][i] = velocity_h.xyz[dim][i] / 2.;
        velocity_h.xyz[dim][i] += force.xyz[dim][i] * dt_md / rmass;   //at time t+0.5dt_md
        velocity.xyz[dim][i] += velocity_h.xyz[dim][i] / 2.;      //at time t
        previous.xyz[dim][i] = pos.xyz[dim][i];               //at time t
        pos.xyz[dim][i] += velocity_h.xyz[dim][i] * dt_md;            //at time t+dt_md
      }
#endif
  }
}


// put the particles back into the box.
void toBox() {

  //change the size of the box:
  for (int dim = 0; dim < 3; dim++)
    xyzL[dim] += xyzL[dim] * vgrad[dim][dim] * dt_md;
  strain_md = (xyzL[0] - xyzL0[0]) / xyzL0[0];

#pragma loop count min(MIN_NP_MD)
#pragma ivdep
  for (int i = 0; i < np_md; i++) {
    //periodicity in z-direction. Change coordinates but not velocities.
    if (pos.xyz[2][i] > 0.5*xyzL[2]) {
      double nz = pos.xyz[2][i] / xyzL[2];
      //     fprintf(stdout,"z-down:%f",pos.xyz[2][i]);
      pos.xyz[2][i] -= ceil(nz) * xyzL[2];
      //fprintf(stdout," %f\n",pos.xyz[2][i]);
    }
    else if (pos.xyz[2][i] < -0.5*xyzL[2]) {
      double nz = - pos.xyz[2][i] / xyzL[2];
      //fprintf(stdout,"z-up:%f",pos.xyz[2][i]);
      pos.xyz[2][i] += ceil(nz) * xyzL[2];
      //fprintf(stdout," %f\n",pos.xyz[2][i]);
    }
    //periodicity in y-direction. Change coordinates but not velocities.
    if (pos.xyz[1][i] > 0.5*xyzL[1]) {
      double ny = pos.xyz[1][i] / xyzL[1];
      //fprintf(stdout,"y-down:%f",pos.xyz[1][i]);
      pos.xyz[1][i] -= ceil(ny) * xyzL[1];
      //fprintf(stdout," %f\n",pos.xyz[1][i]);
    }
    else if (pos.xyz[1][i] < -0.5*xyzL[1]) {
      double ny = - pos.xyz[1][i] / xyzL[1];
      //fprintf(stdout,"y-up:%f",pos.xyz[1][i]);
      pos.xyz[1][i] += ceil(ny) * xyzL[1];
      //fprintf(stdout," %f\n",pos.xyz[1][i]);
    }
    //periodicity in x-direction. Change both coordinates and velocities.
    if (pos.xyz[0][i] > 0.5*xyzL[0]) {
      double nx = pos.xyz[0][i] / xyzL[0];
      //fprintf(stdout,"x-down:%f",pos.xyz[0][i]);
      pos.xyz[0][i] -= ceil(nx) * xyzL[0];
      //fprintf(stdout," %f\n",pos.xyz[0][i]);
      //fprintf(stdout,"x-down vel:%f",velocity.xyz[0][i]);
      velocity.xyz[0][i] -= ceil(nx) * vgrad[0][0];
      //fprintf(stdout," %f\n",velocity.xyz[0][i]);
    }
    else if (pos.xyz[0][i] < -0.5*xyzL[0]) {
      double nx = - pos.xyz[0][i] / xyzL[0];
      //fprintf(stdout,"x-up:%f",pos.xyz[0][i]);
      pos.xyz[0][i] += ceil(nx) * xyzL[0];
      //fprintf(stdout," %f\n",pos.xyz[0][i]);
      //fprintf(stdout,"x-up vel:%f",velocity.xyz[0][i]);
      velocity.xyz[0][i] += ceil(nx) * vgrad[0][0];
      //fprintf(stdout," %f\n",velocity.xyz[0][i]);
    }
  }
}

// this method calculates the virial stresses
void parStress() {

  int nghbr;
  double newforce[3];
  double avgv[3];

#ifdef __INTEL_COMPILER
  stress[:][:] = 0.0;
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      stress[i][j] = 0.0;
#endif

#ifdef __INTEL_COMPILER
# pragma loop count min(MIN_NP_MD)
  avgv[0] = __sec_reduce_add(velocity.xyz[0][0:np_md]) / np_md;
# pragma loop count min(MIN_NP_MD)
  avgv[1] = __sec_reduce_add(velocity.xyz[1][0:np_md]) / np_md;
# pragma loop count min(MIN_NP_MD)
  avgv[2] = __sec_reduce_add(velocity.xyz[2][0:np_md]) / np_md;
#else
  avgv[0] = 0.0;
  avgv[1] = 0.0;
  avgv[2] = 0.0;
# pragma loop count min(MIN_NP_MD)
  for (int ip = 0; ip < np_md; ip++) {
    // particle velocity
    avgv[0] += velocity.xyz[0][ip];
    avgv[1] += velocity.xyz[1][ip];       //vector avgv += vector velocities
    avgv[2] += velocity.xyz[2][ip];
  }
  avgv[0] /= np_md;
  avgv[1] /= np_md;       //vector avgv /= scalar np_md
  avgv[2] /= np_md;
#endif

#pragma loop count min(MIN_NP_MD)
  for (int i = 0; i < np_md; i++) {
    //matrix stress -= scalar rmass * vector (column) (velocities - avgv) * vector (row) (velocities - avgv)
#ifdef __INTEL_COMPILER
    stress[0][:] -= rmass * (velocity.xyz[0][i] - avgv[0]) * (velocity.xyz[:][i] - avgv[:]);
    stress[1][:] -= rmass * (velocity.xyz[1][i] - avgv[1]) * (velocity.xyz[:][i] - avgv[:]);
    stress[2][:] -= rmass * (velocity.xyz[2][i] - avgv[2]) * (velocity.xyz[:][i] - avgv[:]);
#else
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        stress[r][c] -= rmass * (velocity.xyz[r][i] - avgv[r]) * (velocity.xyz[c][i] - avgv[c]);
#endif
    //CHANGE HERE
  }
  double xyzL_local[3] = {xyzL[0], xyzL[1], xyzL[2]};
  printf("kinetic contribution: %f eV\n",stress[0][0]);
  for (int i = 0; i < np_md; i++) {
    double force_prod_i[3][np_md];
    for (int j = 0; j < np_md; j++) {
      if (i == j) continue;
      // normal force
      double dis[3] = {
        (pos.xyz[0][i] - pos.xyz[0][j]),
        (pos.xyz[1][i] - pos.xyz[1][j]),
        (pos.xyz[2][i] - pos.xyz[2][j])
      };
      if (fabs(dis[0]) > 0.5*xyzL_local[0]) {
        double s = 1.0;
        if (dis[0] < 0)
          s = -1.0;
        dis[0] -= s*xyzL_local[0];
      }
      if (fabs(dis[1]) > 0.5*xyzL_local[1]) {
        double s = 1.0;
        if (dis[1] < 0)
          s = -1.0;
        dis[1] -= s*xyzL_local[1];
      }
      if (fabs(dis[2]) > 0.5*xyzL_local[2]) {
        double s = 1.0;
        if (dis[2] < 0)
          s = -1.0;
        dis[2] -= s*xyzL_local[2];
      }

      float distance2 = dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2];
      float inv_distance2 = 1.0f/distance2;
      double sigma6 = pow(sigma, 6.0);
      double sigma12 = sigma6*sigma6;
      double inv_distance8 = pow(inv_distance2, 4);
      double inv_distance14 = pow(inv_distance2, 7);
      double fnorm = 24. * epsilon * (2.0*sigma12 * inv_distance14 - sigma6 * inv_distance8);

#ifdef __INTEL_COMPILER
      newforce[:] = fnorm * dis[:];
      stress[0][:] += 0.5 * (previous.xyz[0][j] - previous.xyz[0][i]) * newforce[:];
      stress[1][:] += 0.5 * (previous.xyz[1][j] - previous.xyz[1][i]) * newforce[:];
      stress[2][:] += 0.5 * (previous.xyz[2][j] - previous.xyz[2][i]) * newforce[:];
#else
      newforce[0] = fnorm * dis[0];
      newforce[1] = fnorm * dis[1];
      newforce[2] = fnorm * dis[2];
      //      if (distance2 > pow(2.55,2) && distance2 < pow(3.6,2)) printf("distance = %f, force = %f\n", sqrt(distance2),sqrt(pow(newforce[0],2)+pow(newforce[1],2)+pow(newforce[2],2)));
      for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
          stress[r][c] += 0.5 * (previous.xyz[r][j] - previous.xyz[r][i]) * newforce[c];
#endif
      
    }
 }

  printf("potential contribution: %f eV\n",stress[0][0]);
#ifdef __INTEL_COMPILER
  stress[:][:] /= __sec_reduce_mul(xyzL[:]);
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      stress[i][j] /= xyzL[0]*xyzL[1]*xyzL[2]; //matrix stress /= scalar xyzL[0]*xyzL[1]*xyzL[2]
#endif

 //  for (int i = 0; i < np_md; i++) {
//     //END CHANGE
//     for (int j = 0; j < i; j++) {
//       double dis[3] = {
//         (pos.xyz[0][i] - pos.xyz[0][j]),
//         (pos.xyz[1][i] - pos.xyz[1][j]),
//         (pos.xyz[2][i] - pos.xyz[2][j])
//       };
//       if(abs(dis[0]) > 0.5*xyzL[0]) dis[0] = -copysign(xyzL[0],dis[0]) + dis[0];   //periodicity
//       if(abs(dis[1]) > 0.5*xyzL[1]) dis[1] = -copysign(xyzL[1],dis[1]) + dis[1];
//       if(abs(dis[2]) > 0.5*xyzL[2]) dis[2] = -copysign(xyzL[2],dis[2]) + dis[2];

//       double distance = sqrt(dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]);
//       double fnorm = 24. * epsilon * (2.*pow(sigma,12) / pow(distance,13)
//                                        - pow(sigma,6) / pow(distance,7));

//       double xyzn[3] = {
//         dis[0]/distance,
//         dis[1]/distance,
//         dis[2]/distance
//       };

//       //matrix stress += 0.5 * vector (column) (previous[i] - previous [j]) * vector (row) newforce
// #ifdef __INTEL_COMPILER
//       newforce[:] = fnorm * xyzn[:];
//       stress[0][:] += 0.5 * (previous.xyz[0][j] - previous.xyz[0][i]) * newforce[:];
//       stress[1][:] += 0.5 * (previous.xyz[1][j] - previous.xyz[1][i]) * newforce[:];
//       stress[2][:] += 0.5 * (previous.xyz[2][j] - previous.xyz[2][i]) * newforce[:];
// #else
//       newforce[0] = fnorm * xyzn[0];
//       newforce[1] = fnorm * xyzn[1];
//       newforce[2] = fnorm * xyzn[2];
//       for (int r = 0; r < 3; r++)
//         for (int c = 0; c < 3; c++)
//           stress[r][c] += 0.5 * (previous.xyz[r][j] - previous.xyz[r][i]) * newforce[c];
// #endif
//     }
//   }
//     printf("potential contribution: %f eV\n",stress[0][0]);
// #ifdef __INTEL_COMPILER
//   stress[:][:] /= __sec_reduce_mul(xyzL[:]);
// #else
//   for (int i = 0; i < 3; i++)
//     for (int j = 0; j < 3; j++)
//       stress[i][j] /= xyzL[0]*xyzL[1]*xyzL[2]; //matrix stress /= scalar xyzL[0]*xyzL[1]*xyzL[2]
// #endif
}

// restart file at the end of relaxation
void dataOutput() {

  FILE *varout;

  varout = fopen("restart.dat","w");
  fprintf(varout,"%d \n",np_md);
  fprintf(varout,"%f %f %f \n",xyzL[0],xyzL[1],xyzL[2]);
  fprintf(varout,"%f \n",rmass);
  for (int i = 0; i < np_md; i++)
    fprintf(varout,"%d %f %f %f\n", i+1, previous.xyz[0][i],previous.xyz[1][i],previous.xyz[2][i]);
  for (int i = 0; i < np_md; i++)
    fprintf(varout,"%d %f %f %f\n", i+1, velocity.xyz[0][i],velocity.xyz[1][i],velocity.xyz[2][i]);
  fclose(varout);

}


void md(int relax, int particle, double SR, double strain) {

  gamxx0 = SR; 
  relaxation = relax;
  final_strain = strain;
  if (relaxation == 1)
    printf("relaxation of MD\n");
  else if (relaxation == 0)
    printf("MD compression\n");
  else {
    printf("unknown value for relax parameter in md\n");
    abort();
  }
    
  par = particle;

  //fprintf(stdout,"setDefault\n");
  setDefault();

  //fprintf(stdout,"dataInput\n");
  dataInput();

  //fprintf(stdout,"setInit\n");
  setInit();
  gamxx0 *= xyzL0[0] * 1e-12; //strain rate in 1/s converted to velocity in A/ps
  time_md = 0.0;

  // printf("calculate nstep\n");
  if (relaxation ==1 )
    nstep = (int) ( (trelax - time_md) / dt_md) + 1;
  else {
    if (abs(SR) <= 1e-3 && abs(final_strain) <= 1e-3)
      nstep = (int) ( (maxtime - time_md) / dt_md) +1;
    else if (abs(SR) <= 1e-3 && abs(final_strain) >= 0.1) {
      printf("too many timesteps, time to compression = %f\n",final_strain / SR);
      abort();
    }
    else
      nstep = (int) ( (maxtime + final_strain/ SR - time_md) / dt_md) +1;
  }
  printf("nstep = %d\n",nstep);

  energy();
  printf("ke = %f; pe = %f; te = %f\n",kenergy,penergy,penergy+kenergy);

  for (it = 0; it < nstep; it++) {
    //fprintf(stdout,"newVgrad\n");
    newVgrad();
    if ((relaxation == 0) && strain_md < final_strain) {
      //fprintf(stdout,"velAdjust\n");
      velAdjust();
    }
    if (time_md == 0 || it % freq == 0) {
      fprintf(stdout,"     Istep = %d, Time = %f\n",it,time_md);
      // FILE *outenergy = fopen("energy.txt","a");
      // if (time_md > 0)
      //   fprintf(outenergy,"%f",time_md - dt_md);
      // else
      //   fprintf(outenergy,"%f",time_md);
      // //fprintf(stdout,"energy\n");
      // energy();

      //fprintf(outenergy," %f %f %f\n",kenergy,penergy,penergy+kenergy);
      //fclose(outenergy);

      //fprintf(stdout,"parStress\n");
      //parStress();

      //fprintf(stdout,"dataOutput\n");
      //dataOutput();
    }
    //fprintf(stdout,"parMove\n");
    parMove();

    //fprintf(stdout,"toBox\n");
    toBox();

    time_md = (it + 1) * dt_md;
  }

  if (relaxation == 1) {
    dataOutput();
  }
  //printf("parStress\n");
  energy();
  printf("ke = %f; pe = %f; te = %f\n",kenergy,penergy,penergy+kenergy);
  parStress();
  //printf("transfering sigma_11\n");
  //sigma_11[par] = 1.60219 * 1e-9 * stress[0][0] * xyzL[1] * xyzL[2];  //1-D stresses in Pa * m^2
  sigma_11[par] = 1.60219 * 1e11 * stress[0][0]; //stresses in Pa
  printf("end of MD, sigma = %lf GPa %lf eV/A^3\n",sigma_11[par]/1e9,stress[0][0]);
}




//-----------------
//MPM
//-----------------
double shape(int inode, double xp) {
  double shape;
  double xpv = xp - nodex[inode];
  if (xpv <= -dx_mpm)
    shape = 0.0;
  else if (xpv <= 0)
    shape = 1.0 + xpv/dx_mpm;
  else if (xpv <= dx_mpm)
    shape = 1.0 - xpv/dx_mpm;
  else
    shape = 0;
  return shape;
}

double derShape(int inode, double xp) {
  // flip part
  double der;
  double xi = nodex[inode];
  double xm1 =-dx_mpm;
  double xp1 = length_mpm+dx_mpm;

  if(inode>=1)
    xm1 = nodex[inode-1];

  if(inode<=nnodes-2)
    xp1 = nodex[inode+1];
  if(fabs(xp-xi)> 2.0*dx_mpm ) return 0.0;
        
  //original MPM part
  double der1;
  if (xp <= xm1) der1 = 0.0;
  else if (xp < xi) der1 = 1.0/dx_mpm;
  else if (xp < xp1) der1 = -1.0/dx_mpm;
  else  der1 = 0.0;

  return der1;
}

double *parToGridForce() {
  double *force = new double[nnodes];
  for (int ip = 0; ip < np_mpm; ip++) {
    for (int inode = 0; inode < nnodes; inode++) {
      force[inode] -= 2.0*lp[ip] * sigma_11[ip]
	* derShape(inode, xp[ip]);
    }
  }
  return force;
}

void parToGridMass() {
  double sum = 0.0;
  for (int inode = 0; inode < nnodes; inode++) {
    sum = 0.0;
    for (int ip = 0; ip < np_mpm; ip++)
      sum += shape(inode, xp[ip]) * mp_mpm[ip];
    nodeMass[inode] = sum;
  }
}

void addParticleAcceleration() {
  double *parForce = parToGridForce();
  for (int inode = 0; inode < nnodes; inode++) {
    langVel[inode] += parForce[inode] / (nodeMass[inode] + 1e-8) * dt_mpm;
  }
}

void addBodyForce() {
  gravity = targetG * time_mpm /loadingTime;
  if (gravity < targetG)
    gravity = targetG;
  for (int inode = 0; inode < nnodes; inode++){
    if(nodeMass[inode] > 0.0)	langVel[inode] += gravity * dt_mpm;
  }	
  //boundary conditions
  langVel[0] = bcfirst;
  langVel[nnodes-1] = bclast;
}

void gridToParStress() {
  double strainRate;		
  for (int ip = np_mpm-1; ip >= 0 ; ip--) {
    strainRate = 0.0;
    for (int inode = 0; inode < nnodes; inode++) {
      strainRate += langVel[inode] * derShape(inode, xp[ip]);
    }
    //strain calculation
    strain[ip] += strainRate * dt_mpm;
    printf("calling MD\n");
    printf("ip = %d strainRate = %f strain = %f\n",ip,strainRate,strain[ip]);
    //call MD here on each particle
    md(0,ip,strainRate,strain[ip]);
    //sigma_11[ip] += YoungM * strainRate * dt_mpm;
    lp[ip] *= (1.0 + strainRate * dt_mpm);
  }
}

void gridToParVel() {
  for (int ip = 0; ip < np_mpm; ip++) {
    for (int inode = 0; inode < nnodes; inode++) {
      velp[ip] += (langVel[inode] - nodeVel[inode])
	* shape(inode, xp[ip]);
    }
  }
}

void movePars() {
  for (int ip = 0; ip < np_mpm; ip++) {
    for (int inode = 0; inode < nnodes; inode++) {
      xp[ip] += 0.5 * (langVel[inode] + nodeVel[inode])
	* shape(inode, xp[ip]) * dt_mpm;
    }
  }
}

void parToGridvel() {

  double sum = 0.0;
  for (int inode = 0; inode < nnodes; inode++) {
    sum = 0.0;
    for (int ip = 0; ip < np_mpm; ip++)
      sum += shape(inode, xp[ip]) * mp_mpm[ip] * velp[ip];
    langVel[inode]=
      nodeVel[inode] = sum / (nodeMass[inode] + 1.e-8);			
  }
  //boundary conditions
  langVel[nnodes-1] = bclast;
  nodeVel[nnodes-1] = bclast;
  langVel[0] = bcfirst;
  nodeVel[0] = bcfirst;
}

void outputData(int num) {
  char filename[20],filename2[20],str1[6];

  sprintf(str1,"%d",num);
  strcpy(filename,"MPM_");
  strncat(filename,str1,6);//change if num is more than 6 digits
  strcat(filename,".dat");
  sprintf(str1,"%d",num);
  strcpy(filename2,"MPM2_");
  strncat(filename2,str1,6);//change if num is more than 6 digits
  strcat(filename2,".dat");
  FILE *plotfile = fopen(filename,"w");
  for (int ip = 0; ip < np_mpm; ip++) {
    fprintf(plotfile,"%f  %f  %f\n",xp[ip],sigma_11[ip],velp[ip]);
  }
  FILE *plotfile2 = fopen(filename2,"w");
  for (int ip = 0; ip < nnodes; ip++) {
    fprintf(plotfile2,"%f  %f  %f\n",nodex[ip],nodeVel[ip],langVel[ip]);
  }
  fclose(plotfile);
  fclose(plotfile2);
  //error handling needed here
}

void input() {
  timeOut = 1.0e-4;//0.0012;
  dt_mpm = 1.0e-4; // sec
  length_mpm = 0.5; // m
  ncells = 50;
  nnodes = ncells + 1;
  double rho = 8.94 * 1e3; // kg/m^3
  dx_mpm = length_mpm / ncells;  //m

  //YoungM = 1.0e6;

  parsPerCell = 2;

  parSize = dx_mpm / parsPerCell;  //m

  parVol = dx_mpm / parsPerCell;  //m

  nodex = new double[nnodes];
  nodeVel = new double[nnodes];
  langVel = new double[nnodes];
  nodeMass = new double[nnodes];
  grdVol = new double[nnodes];

  for (int inode = 0; inode < nnodes; inode++) {
    nodex[inode] = dx_mpm * inode;  //m
    nodeVel[inode] = 0.0;    
    grdVol[inode] = dx_mpm;     //m
  }
  
  np_mpm = (ncells) * parsPerCell;
  xp = new double[np_mpm];
  l0_mpm= new double[np_mpm];
  lp = new double[np_mpm];
  sigma_11 = new double[np_mpm];
  strain = new double[np_mpm];
  mp_mpm = new double[np_mpm];
  velp = new double[np_mpm];
  for (int ip = 0; ip < np_mpm; ip++) {
    xp[ip] = dx_mpm / parsPerCell / 2.0 + dx_mpm / parsPerCell * ip;    //m
    mp_mpm[ip] = rho * dx_mpm / parsPerCell;  //kg
    lp[ip] = parSize / 2.0;    //m
    sigma_11[ip] = 0;
    strain[ip] = 0;
  }
  //initial conditions
  nodeVel[nnodes-1] = iclast;
  langVel[nnodes-1] = iclast;
  velp[np_mpm-1] = 0.75*iclast;
  velp[np_mpm-2] = 0.25*iclast;
  nodeVel[0] = icfirst;
  langVel[0] = icfirst;
  velp[0] = 0.75*icfirst;
  velp[1] = 0.25*icfirst;

  parToGridMass();
}

void mainLoop() {
  int num = 1;
  outputData(num++);
  //printf("Relaxing MD\n");
  md(1,0,0,0); //relaxation of MD model
  for (int ip = 1; ip < np_mpm; ip++) {
    sigma_11[ip] = sigma_11[0];
  }
  while (time_mpm < loadingTime) {
    //printf("addParticleAcceleration\n");
    addParticleAcceleration();
    //printf("addBodyForce\n");
    addBodyForce();
    //printf("gridToParStress\n");
    gridToParStress();
    //printf("gridToParVel\n");
    gridToParVel();
    //printf("movePars\n");
    movePars();
    //printf("parToGridMass\n");
    parToGridMass();
    //printf("parToGridvel\n");
    parToGridvel();

    if (isnan(xp[0])) {
      printf("NaN at time = %f\n",time_mpm);
      exit(1);
    }
    time_mpm += dt_mpm;
            
    sinceLastOut += dt_mpm;
			
    if(fabs(sinceLastOut - timeOut) < dt_mpm)
      {
	sinceLastOut = 0.0;
	outputData(num++);
	printf( "time_mpm = %f\n",time_mpm);
      }
  }

  printf("The End\n");
}

int main() {
  input();
  mainLoop();  
}
