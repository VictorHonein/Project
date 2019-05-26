#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/swirl.h"
#include "two-phase.h"
#include "tension.h"
#include "axistream.h"
#include "view.h"

u.t[left] = dirichlet(0);
w[left]   = dirichlet(y);

u.t[right] = dirichlet(0);
w[right]   = dirichlet(0);

u.t[top] = dirichlet(0);
w[top]   = dirichlet(0);

double G = 0.25,
       Fr = 0.88,
       Re = 3063,
       We = 3153,
       rhor = 1.205/1.2107e3,
       mur = 18.2e-6/6.09e-2;

int main()
{
  N = 128;
  mu1 = 1./Re;
  rho2 = rho1*rhor;
  mu2 = mu1*mur;
  //  f.sigma = 1./We;
  DT = 0.1;
  run();
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1./Fr;
}

event init (i = 0) {
  fraction (f, G - x);
}

#if 0
event snapshots (i += 100) {
  scalar psi[];
  psi[left] = dirichlet (0);
  psi[right] = dirichlet (0);
  psi[top] = dirichlet (0);
  axistream (u, psi);
  dump();
}
#endif

event logfile (i += 10)
{
  double ke = 0.;
  foreach (reduction(+:ke))
    ke += dv()*rho(f[])*(sq(u.x[]) + sq(u.y[]) + sq(w[]));
  fprintf (stderr, "%g %g %g %g\n", t, dt, ke, statsf (f).sum);
}

event movie (t += 1; t <= 300)
{
  view (fov = 22.7232, quat = {0,0,-0.707,0.707},
	tx = -0.12824, ty = -0.446981,
	width = 830, height = 432);
  box();

scalar psi[];
  psi[left] = dirichlet (0);
  psi[right] = dirichlet (0);
  psi[top] = dirichlet (0);
  axistream (u, psi);

 squares ("psi", linear = true, spread = -1);
  isoline ("psi", n = 21, min = 0, max = 0.006);
  isoline ("psi", -0.0012, lc = {0,1,0});
  draw_vof ("f", lc = {1,0,0}, lw = 2);
  translate (y = -1.15) {
    box();
    squares ("w", linear = true, spread = -1);
    isoline ("w", n = 21, min = 0, max = 1);
    draw_vof ("f", lc = {1,0,0}, lw = 2);
  }
  save ("movie.mp4");
}

#if 0 // TREE
event adapt (i++) {
  adapt_wavelet ({u,w}, (double[]){1e-3,1e-3,1e-3}, minlevel = 5, maxlevel = 8);
}
#endif
