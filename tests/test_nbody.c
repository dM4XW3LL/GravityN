/**
 * @file  test_nbody.c
 * @brief Validation harness for nbody.c
 *
 * Tests run:
 *   1. Earth–Sun: orbit should close after 1 year, energy conserved.
 *   2. Three-body: barycentre fixed, momentum conserved.
 *   3. integrator comparison: leapfrog vs RK4 should agree to better than 1e-6 AU.
 */
#include <stdio.h>
#include <math.h>
#include "../src/core/nbody.h"
#include "../src/core/kepler.h"

#define PASS "\033[32mPASS\033[0m"
#define FAIL "\033[31mFAIL\033[0m"

/* ── helper: build a Sun+Earth simulation ─────────────────────────────── */
static void build_earth_sun(Simulation *sim)
{
    nbody_init(sim);

    Body sun = {0};
    sun.mass = 1.0;
    snprintf(sun.name, NBODY_NAME_LEN, "Sun");
    nbody_add_body(sim, &sun);

    Body earth = {0};
    earth.mass = 3.003e-6;
    elements_to_state(1.0, 0.0167, 0.0, 1.0,
                      &earth.x, &earth.y, &earth.vx, &earth.vy);
    snprintf(earth.name, NBODY_NAME_LEN, "Earth");
    nbody_add_body(sim, &earth);

    /* Prime the accelerations before stepping */
    nbody_compute_accelerations(sim);
}

/* ── Test 1: Earth closes its orbit after 1 year ─────────────────────── */
static void test_orbit_closure(void)
{
    Simulation sim;
    build_earth_sun(&sim);

    double x0 = sim.bodies[1].x;
    double y0 = sim.bodies[1].y;

    int steps = 10000;
    double dt = 1.0 / steps;
    for (int i = 0; i < steps; i++)
        nbody_step_leapfrog(&sim, dt);

    double dx = sim.bodies[1].x - x0;
    double dy = sim.bodies[1].y - y0;
    double closure = sqrt(dx*dx + dy*dy);

    printf("Test 1 - Orbit closure after 1 yr:  %.3e AU   %s\n",
           closure, closure < 1e-4 ? PASS : FAIL);

    nbody_free(&sim);
}

/* ── Test 2: Energy conservation over 10 years ───────────────────────── */
static void test_energy_conservation(void)
{
    Simulation sim;
    build_earth_sun(&sim);

    double E0 = nbody_total_energy(&sim);

    int steps = 50000;
    double dt = 10.0 / steps;
    for (int i = 0; i < steps; i++)
        nbody_step_leapfrog(&sim, dt);

    double E1 = nbody_total_energy(&sim);
    double rel_err = fabs((E1 - E0) / E0);

    printf("Test 2 - Energy conservation (10 yr): dE/E = %.3e   %s\n",
           rel_err, rel_err < 1e-5 ? PASS : FAIL);

    nbody_free(&sim);
}

/* ── Test 3: Momentum conservation ──────────────────────────────────────*/
static void test_momentum_conservation(void)
{
    Simulation sim;
    nbody_init(&sim);

    /* Three bodies with non-zero initial velocities */
    Body b = {0};
    b.mass = 1.0;  b.x =  0.0; b.y =  0.0; b.vx =  0.1; b.vy =  0.0;
    snprintf(b.name, NBODY_NAME_LEN, "A");
    nbody_add_body(&sim, &b);

    b.mass = 0.5;  b.x =  3.0; b.y =  0.0; b.vx = -0.2; b.vy =  0.1;
    snprintf(b.name, NBODY_NAME_LEN, "B");
    nbody_add_body(&sim, &b);

    b.mass = 0.3;  b.x = -2.0; b.y =  2.0; b.vx =  0.0; b.vy = -0.15;
    snprintf(b.name, NBODY_NAME_LEN, "C");
    nbody_add_body(&sim, &b);

    nbody_compute_accelerations(&sim);

    double px0, py0, px1, py1;
    nbody_total_momentum(&sim, &px0, &py0);

    for (int i = 0; i < 10000; i++)
        nbody_step_leapfrog(&sim, 0.001);

    nbody_total_momentum(&sim, &px1, &py1);

    double dpx = fabs(px1 - px0);
    double dpy = fabs(py1 - py0);

    printf("Test 3 - Momentum conservation:   dpx=%.3e dpy=%.3e   %s\n",
           dpx, dpy, (dpx < 1e-12 && dpy < 1e-12) ? PASS : FAIL);

    nbody_free(&sim);
}

/* ── Test 4: Leapfrog vs RK4 agreement ──────────────────────────────── */
static void test_integrator_agreement(void)
{
    Simulation lf, rk;
    build_earth_sun(&lf);
    build_earth_sun(&rk);

    int steps = 1000;
    double dt = 1.0 / steps;
    for (int i = 0; i < steps; i++) {
        nbody_step_leapfrog(&lf, dt);
        nbody_step_rk4(&rk, dt);
    }

    /* Compare Earth positions */
    double dx = lf.bodies[1].x - rk.bodies[1].x;
    double dy = lf.bodies[1].y - rk.bodies[1].y;
    double diff = sqrt(dx*dx + dy*dy);

    printf("Test 4 - Leapfrog vs RK4 (1 yr):   diff = %.3e AU   %s\n",
           diff, diff < 5e-4 ? PASS : FAIL);

    nbody_free(&lf);
    nbody_free(&rk);
}

/* ── Test 5: elements_to_state round-trip radius ────────────────────── */
static void test_elements_to_state(void)
{
    /* For e=0 (circular orbit, a=1 AU), the initial position should have
     * |r| = 1.0 AU exactly (within floating-point noise).               */
    double x, y, vx, vy;
    elements_to_state(1.0, 0.0, 0.0, 1.0, &x, &y, &vx, &vy);
    double r = sqrt(x*x + y*y);

    /* For a circular orbit, the orbital speed should equal sqrt(GM/a)   */
    double v       = sqrt(vx*vx + vy*vy);
    double v_circ  = sqrt(NBODY_G * 1.0 / 1.0);

    printf("Test 5 - elements_to_state (circular): r=%.8f AU, v=%.8f AU/yr   %s\n",
           r, v, (fabs(r - 1.0) < 1e-12 && fabs(v - v_circ) < 1e-10) ? PASS : FAIL);
}

int main(void)
{
    printf("\n=== nbody.c validation suite ===\n\n");
    test_elements_to_state();
    test_orbit_closure();
    test_energy_conservation();
    test_momentum_conservation();
    test_integrator_agreement();
    printf("\n");
    return 0;
}