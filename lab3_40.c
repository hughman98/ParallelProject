/* Name: J. Hugh Wright
 * Date: 3/19/19
 * Description: This program simulates a universe of 10,000 bodies spread out across 40 processors.
 * 		It uses upc for parallelization and was created as a learning exercise for a college course.
 */

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "upc_relaxed.h"


#define N THREADS*250
#define G 6.67e-11
#define TIMESTEP 0.25
#define NSTEPS 10

/* 
 * body data structure
 */
struct body_s {
	double x;
	double y;
	double z;
	double dx;
	double dy;
	double dz;
	double mass;
};
typedef struct body_s body_t;


/* 
 * function prototypes
 */
void init(void);
double dist(double dx, double dy, double dz);
 
shared int globalarray[THREADS];
shared [250] body_t bodies[N];  // array of N-bodies at timestep t
shared [250] body_t next[N];    // array of N-bodies at timestep t+1

int eprintf(const char *format, ...) {
  va_list ap;
  int ret;

  if (MYTHREAD == 0) {
    va_start(ap, format);
    ret = vfprintf(stdout, format, ap);
    va_end(ap);
    return ret;
  }
  else
    return 0;
}

/**
 * init - give the planets initial values for position, velocity, mass
 */
void init(void) {
	upc_forall (int i=0; i<N; i++; &bodies[i]) {
		bodies[i].x = 100.0 * (i + 0.1);
		bodies[i].y = 200.0 * (i + 0.1);
		bodies[i].z = 300.0 * (i + 0.1);
		bodies[i].dx = i + 400.0;
		bodies[i].dy = i + 500.0;
		bodies[i].dz = i + 600.0;
		bodies[i].mass = 10e6 * (i + 100.2);
	}
}



/**
 * dist - determine the distance between two bodies
 *    @param dx - distance in the x dimension
 *    @param dy - distance in the y dimension
 *    @param dz - distance in the z dimension
 *    @return distance 
 */
double dist(double dx, double dy, double dz) {
	return sqrt((dx*dx) + (dy*dy) + (dz*dz));;
}



/**
 *  print_body - prints a body for debugging
 *    @param b - body to print
 */
void print_body(body_t b) {
	printf("x: %7.3f y: %7.3f z: %7.3f dx: %7.3f dy: %7.3f dz: %7.3f\n",
			b.x, b.y, b.z, b.dx, b.dy, b.dz);
}

/*
 * get_wctime - returns wall clock time as double
 *   @return double representation of wall clock time
 */
double get_wctime(void) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (tv.tv_sec + 1E-6 * tv.tv_usec);
}



/**
 * main
 */

int main(int argc, char **argv) {
  setenv("PTL_IGNORE_UMMUNOTIFY","1",1);
	double start, tsstart;

  setbuf(stdout, NULL);

	// setup initial conditions
  init();
  /*
  //Print out threads in order.
  for (int i=0; i<THREADS; i++) {
    if (i == MYTHREAD){
      char buf[100];
      gethostname(buf, sizeof(buf));
      printf("UPC thread %d of %d on %s\n",
          MYTHREAD, THREADS, buf);
      globalarray[i] = i;
    }
    upc_barrier;
    //sleep(.01);
  } */

  if (MYTHREAD == 0)
	  printf("beginning N-body simulation of %d bodies with %d processes.\n", N, THREADS);

  upc_barrier;
  //sleep(1); //to make sure the prints go in the right order

	start = get_wctime();

  double d, f;       // distance, force
  double dx, dy, dz; // position deltas
  double fx, fy, fz; // force components
  double ax, ay, az; // acceleration components

  // for each timestep in the simulation
  for (int ts=0; ts<NSTEPS; ts++) {

    tsstart = get_wctime();

    body_t l_bodies[250]; // Change this to = 10000/THREADS for the final run
      
    // Computes the force for all the objects in the universe.
    upc_forall (int i=0; i<N; i++; &bodies[i]) {
  
      fx = fy = fz = 0.0;     
      
      for (int j=0; j<THREADS; j++) {
        
        upc_memget( &l_bodies[0], &bodies[j*250], sizeof(body_t)*250 ); // CHANGE THIS WHEN EDITING THREAD SIZE

        for(int k=0; k<250; k++){ //CHANGE THIS WHEN EDITING THREAD SIZE

          dx = bodies[i].x - l_bodies[k].x;
          dy = bodies[i].y - l_bodies[k].y;
          dz = bodies[i].z - l_bodies[k].z;

          d = dist(dx, dy, dz);

          if (d != 0) {
            f = (G * bodies[i].mass * l_bodies[k].mass) / (d*d);

            fx += (f * dx) /d;
            fy += (f * dy) /d;
            fz += (f * dz) /d;

          }
        }
      }
      
      ax = fx / bodies[i].mass;
      ay = fy / bodies[i].mass;
      az = fz / bodies[i].mass;     

      next[i].dx = bodies[i].dx + (TIMESTEP * ax);
      next[i].dy = bodies[i].dy + (TIMESTEP * ay);
      next[i].dz = bodies[i].dz + (TIMESTEP * az);

      next[i].x = bodies[i].x + (TIMESTEP * bodies[i].dx);
      next[i].y = bodies[i].y + (TIMESTEP * bodies[i].dy);
      next[i].z = bodies[i].z + (TIMESTEP * bodies[i].dz);

      next[i].mass = bodies[i].mass;
    }
    
    // copy the t+1 state to be the new time t
    upc_forall (int i=0; i< N; i++; &bodies[i]) {
      upc_memcpy(&bodies[i], &next[i], sizeof(body_t));
    }

    upc_barrier;
    if (MYTHREAD == 0)
      printf("timestep %d complete: %7.3f ms\n", ts, (get_wctime()-tsstart)*1000);
  }
  if (MYTHREAD == 0){
    printf("simulation complete: %9.3f ms\n", (get_wctime()-start)*1000);
  }
	return 0;
}
