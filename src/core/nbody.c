/**
 * @file nbody.c
 * @author Diogo Pereira Gomes (https://github.com/dM4XW3LL)
 * @brief GravityN - N-body gravitational simulator implementation
 * 
 * @details
 * This file implements every function declared in nbody.h. Read the header
 * for the physical model, unit system, and usage examples; this file focuses on
 * *how* each function works rather than *what* it does.
 * 
 * ### Design Philosophy
 * 
  * - **No dynamic allocation in the hot path.**  All body storage lives in the
 *   @c Simulation struct's embedded array, so the integrator never calls
 *   @c malloc / @c free during a step.  This keeps the inner loop cache-friendly
 *   and avoids allocation failures and memory leaks at runtime.
 *
 * - **Newton's-third-law symmetry.**  The acceleration loop only visits each
 *   pair (i, j) with i < j and adds the contribution to both bodies at once,
 *   halving the number of square-root evaluations — the most expensive
 *   operation in the engine.
 *
 * - **Leapfrog is the default.**  The velocity-Verlet form used here stores
 *   the acceleration at the *current* position inside each Body struct, so
 *   a single force evaluation per step is sufficient.  RK4 requires four
 *   evaluations and temporary state copies, which is why it is only used
 *   when explicitly requested. Besides this Leapfrog is Symplectic, thus
 *   conserving the physical constants and properties of the system, stopping
 *   it from diverging and giving very innacurate results with time.
 * 
 * ### Compilation
 * 
 Link against the math library:
 * @code{.sh}
 *   # Shared library (for Python ctypes)
 *   gcc -O2 -fPIC -shared -o nbody.so nbody.c kepler.c -lm
 *
 *   # Static object (for a standalone C program)
 *   gcc -O2 -c nbody.c -o nbody.o
 * @endcode
 * 
 * @version 1.0
 * @date 2026-06-09
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include <math.h>
#include <string.h>   /* memset, memmove */
#include <stdio.h>    /* (not used at runtime; available for debug prints) */
#include "nbody.h"
#include "kepler.h"   /* kepler_newton, eccentric_to_true — used in elements_to_state */

/* ─────────────────────────────────────────────────────────────────────────────
   Internal helpers (not part of the public API)
   ─────────────────────────────────────────────────────────────────────────────*/

/**
 * @internal
 * @brief  Compute gravitational accelerations for an arbitrary body array.
 *
 * @details
 * This is the *real* acceleration kernel, factored out so that both the
 * public @ref nbody_compute_accelerations (which operates on sim->bodies) and
 * the RK4 helper (which operates on temporary copies) can call it without
 * duplicating code.
 *
 * The pairwise loop visits every unique pair (i, j) with i < j exactly once
 * and exploits Newton's third law:
 *
 * @f[
 *   \mathbf{a}_{ij} = -\mathbf{a}_{ji}
 * @f]
 *
 * so the force between bodies i and j contributes simultaneously to
 * bodies[i].ax and bodies[j].ax (with opposite signs).
 *
 * The softened distance is:
 * @f[
 *   d_s = \left(dx^2 + dy^2 + \varepsilon^2\right)^{3/2}
 * @f]
 *
 * @param[in,out] bodies  Array of bodies whose ax/ay fields are overwritten.
 * @param[in]     n       Number of active bodies in the array.
 */
static void compute_accel_array(Body *bodies, int n)
{

    int i, j;
    double dx, dy, dist2, dist3, f;

    /* Zero all accelerations before accumulating contribution.
     * We must not add to stale values from the previous step. */
    for (i = 0; i < n; i++){
        bodies[i].ax = 0.0;
        bodies[i].ay = 0.0;
    }

    /* Pairwise force loop — O(N²), but we visit each pair once only.   */

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {

            /* Vector from body i to body j.                            */

            dx = bodies[j].x - bodies[i].x;
            dy = bodies[j].y - bodies[i].y;

            /* Softened squared distance, then raise to the 3/2 power. 
             * Using pow(x,1.5) is equivalent to doing x*sqrt(x) and is
             * optimized by most compilers with -O2                     */

            dist2 = dx*dx + dy*dy + NBODY_SOFTENING * NBODY_SOFTENING;
            dist3 = pow(dist2, 1.5);

            /* Acceleration magnitude factor: G/dist3.
             * The body masses are multiplied separately below so that 
             * we can apply the force to both bodies with a single divide.  */
            
            f = NBODY_G/dist3;
            
            /* Body i is accelerated *towards* body j by G*m_j/dist^3.   */
            bodies[i].ax += f*bodies[j].mass*dx;
            bodies[i].ay += f*bodies[j].mass*dy;

            /* Body j is accelerated *towards body i by G*m_i/dist^3    */
            /* Opposite reaction -> Newton's third law                  */
            bodies[j].ax += f*bodies[i].mass*dx;
            bodies[j].ay += f*bodies[i].mass*dy;
        }

    }

}

/* ─────────────────────────────────────────────────────────────────────────────
   Simulation lifecycle
   ─────────────────────────────────────────────────────────────────────────────*/

void nbody_init(Simulation *sim)
{
    /*
     * Zero the entire struct.  This correctly initialises:
     *   - all body positions, velocities, and accelerations to 0.0
     *   - all body masses to 0.0  (caller must set mass > 0 before adding)
     *   - all name strings to NUL
     *   - n to 0
     *   - t to 0.0
     */
    (void)sim; /* supresses "unused parameter" warning*/


}

int nbody_add_body(Simulation *sim, const Body *body)
{

    if (sim->n >= NBODY_MAX_BODIES)
        return -1; /* no room, max bodies utilized*/

    /* Copy the source Body by value into the next free slot.
     * The caller's Body can be a stack variable - we don't hold a pointer
     * to it after this function returns.
     */

    sim->bodies[sim->n] = *body;
    sim->n++;
    return 0;
}

int nbody_remove_body(Simulation *sim, int idx)
{

    if (idx < 0 || idx >= sim->n)
        return -1; /*out of range*/
    
    /*
     * Shift every body above the removed body down by one slot
     * memmove handles the overlapping case correctly (which memcpy
     * does not guarantee)
     */

    int bodies_to_move = sim->n - idx - 1;
    if (bodies_to_move > 0){
        memmove(&sim->bodies[idx],
                &sim->bodies[idx+1],
                (size_t)bodies_to_move*sizeof(Body));
    }

    sim->n--;
    return 0;
}

/* ─────────────────────────────────────────────────────────────────────────────
   Initial-condition helpers
   ─────────────────────────────────────────────────────────────────────────────*/

void elements_to_state(double a, double e, double M0, double M_star, double *x, double *y, double *vx, double *vy)
{

    /*
     * Step 1: Solve Kepler's equation M = E - e*sin(E) for the eccentric
     * anomaly E at the initial mean anomaly M0.
     *
     * We borrow kepler_newton() directly from kepler.c.  The iter_count
     * parameter is a remnant of its original benchmarking role; we pass
     * a local int that we discard.
     */
    int iters = 0;
    double E = kepler_newton(M0, e, &iters);

    /*
     * Step 2: Convert eccentric anomaly to Cartesian position.
     *
     * In the orbital plane with the focus (star) at the origin, the
     * position is:
     *
     *   x = a * (cos(E) - e)
     *   y = a * sqrt(1 - e²) * sin(E)
     *
     * Derivation: The parametric ellipse is centred at the *geometric*
     * centre, offset from the focus by a*e.  Shifting the origin to the
     * focus gives the -e term in x.
     */

    double cos_E = cos(E);
    double sin_E = sin(E);

    *x = a*(cos_E-e);
    *y = a*sqrt(1.0-e*e)*sin_E;

    /*
     * Step 3: Compute the orbital velocity.
     *
     * Differentiating the position equations with respect to time and
     * substituting Ė = n / (1 - e*cos(E))  (where n = 2π/P is the mean
     * motion and P = 2π*sqrt(a³ / G*M_star) by Kepler's third law):
     *
     *   vx = -n*a*sin(E) / (1 - e*cos(E))
     *   vy =  n*a*sqrt(1-e²)*cos(E) / (1 - e*cos(E))
     *
     * We compute n directly from Kepler's third law in solar units:
     * n = sqrt(G*M_star / a³) = 2π * sqrt(M_star / a³)
     *
     * The denominator (1 - e*cos(E)) is the ratio dt/dE and accounts for
     * the fact that the body moves faster near periapsis.
     */

    double n = sqrt(NBODY_G * M_star / (a * a * a));
    
    *vx = -n * a * sin_E / (1.0 - e * cos_E);
    *vy = -n * a * sqrt(1.0 - e * e)*cos_E / (1.0 - e * cos_E);

}

/* ─────────────────────────────────────────────────────────────────────────────
   Core physics
   ─────────────────────────────────────────────────────────────────────────────*/

void nbody_compute_accelerations(Simulation *sim)
{

/*
     * Thin public wrapper around the internal helper.
     * Exposed in the API so that the Python side can force an acceleration
     * recompute (e.g. right after adding a new body mid-simulation) without
     * having to take an integration step.
     */

    compute_accel_array(sim->bodies, sim->n);

}

/* ─────────────────────────────────────────────────────────────────────────────
   Leapfrog integrator
   ─────────────────────────────────────────────────────────────────────────────*/


void nbody_step_leapfrog(Simulation *sim, double dt)
{

    /*
     * Velocity-Verlet form of the leapfrog integrator.
     *
     * The key insight behind leapfrog is that velocities and positions are
     * defined at *interleaved* (staggered) time points.  Written out, one
     * full timestep from t to t+dt looks like:
     *
     *   [1] Half velocity kick:   v(t + dt/2) = v(t)        + a(t) * dt/2
     *   [2] Full position drift:  r(t + dt)   = r(t)        + v(t+dt/2) * dt
     *   [3] Recompute forces:     a(t + dt)   = f(r(t+dt)) / m
     *   [4] Half velocity kick:   v(t + dt)   = v(t+dt/2)  + a(t+dt) * dt/2
     *
     * Why is this symplectic?  Because each of the three operations (kick,
     * drift, kick) is itself a shear transformation in phase space that
     * preserves volume (Liouville's theorem).  Composing volume-preserving
     * maps gives a volume-preserving map — the hallmark of a symplectic
     * integrator.
     *
     * Bootstrapping: on the very first call, ax and ay are 0 (set by
     * nbody_init via memset), so step [1] is a no-op and the integrator
     * starts correctly from rest.  However, it is better practice to call
     * nbody_compute_accelerations once before the first step so that the
     * initial forces reflect the actual body positions.
     */

    int i;
    double half_dt = 0.5*dt;

    /* ── Step 1: half velocity kick ─────────────────────────────────────── */

    for (i = 0; i<sim->n; i++){
        sim->bodies[i].vx += sim->bodies[i].ax * half_dt;
        sim->bodies[i].vy += sim->bodies[i].ay * half_dt;
    }

    /* ── Step 2: full position drift ────────────────────────────────────── */

    for (i = 0; i<sim->n; i++){
        sim->bodies[i].x += sim->bodies[i].vx * dt;
        sim->bodies[i].y += sim->bodies[i].vy * dt;
    }

    /* ── Step 3: recompute accelerations at new positions ──────────────── */
    compute_accel_array(sim->bodies,sim->n);

    /* ── Step 4: second half velocity kick ──────────────────────────────── */
    for (i = 0; i<sim->n; i++){
        sim->bodies[i].vx += sim->bodies[i].ax * half_dt;
        sim->bodies[i].vy += sim->bodies[i].ay * half_dt;
    }

    /* Advance the simulation clock. */

    sim->t += dt;

}


/* ─────────────────────────────────────────────────────────────────────────────
   RK4 integrator
   ─────────────────────────────────────────────────────────────────────────────*/

void nbody_step_rk4(Simulation *sim, double dt)
{

    /*
     * Classical 4th-order Runge–Kutta (RK4) integrator.
     *
     * The state of the system is the concatenation of all positions and
     * velocities: y = (r_0, v_0, r_1, v_1, …, r_{n-1}, v_{n-1}).
     * The derivative dy/dt = (v_0, a_0(r), v_1, a_1(r), …).
     *
     * RK4 evaluates this derivative at four carefully chosen intermediate
     * points and combines them to cancel all error terms up to O(dt⁴):
     *
     *   k1 = f(y)
     *   k2 = f(y + dt/2 * k1)
     *   k3 = f(y + dt/2 * k2)
     *   k4 = f(y + dt   * k3)
     *
     *   y(t+dt) = y(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
     *
     * Implementation strategy:
     * We represent the "state + derivative" as a temporary Body array.
     * Each of the four stages is a snapshot of (position, velocity) from
     * which we compute the derivative (velocity, acceleration).
     *
     * Naming convention used below:
     *   b0   = bodies at the start of the step (read-only reference)
     *   tmp  = temporary state at an intermediate RK stage
     *   k1…k4 = derivative arrays (velocity and acceleration at each stage)
     */

    int i, n;
    n = sim->n;

    /* Snapshot of the starting state — we must not modify sim->bodies
     * until the final weighted update at the very end.                      */

    Body b0[NBODY_MAX_BODIES];
    memcpy(b0, sim->bodies, (size_t)n * sizeof(Body));

    /* Derivative arrays: each k holds (vx, vy) as the position derivative
     * and (ax, ay) as the velocity derivative for every body.               */

    double k1vx[NBODY_MAX_BODIES], k1vy[NBODY_MAX_BODIES];
    double k1ax[NBODY_MAX_BODIES], k1ay[NBODY_MAX_BODIES];

    double k2vx[NBODY_MAX_BODIES], k2vy[NBODY_MAX_BODIES];
    double k2ax[NBODY_MAX_BODIES], k2ay[NBODY_MAX_BODIES];
 
    double k3vx[NBODY_MAX_BODIES], k3vy[NBODY_MAX_BODIES];
    double k3ax[NBODY_MAX_BODIES], k3ay[NBODY_MAX_BODIES];
 
    double k4vx[NBODY_MAX_BODIES], k4vy[NBODY_MAX_BODIES];
    double k4ax[NBODY_MAX_BODIES], k4ay[NBODY_MAX_BODIES];

    Body tmp[NBODY_MAX_BODIES];


    /* ── Stage 1: derivative at the current state ──────────────────────── */
    /*
     * The derivative of position is velocity; the derivative of velocity is
     * acceleration.  We compute accelerations at the starting positions.
     */

    compute_accel_array(sim->bodies, n); /* fills bodies[i].ax/ay */

    for (i = 0; i < n; i++){
        k1vx[i] = sim->bodies[i].vx;
        k1vy[i] = sim->bodies[i].vy;
        k1ax[i] = sim->bodies[i].ax;
        k1ay[i] = sim->bodies[i].ay;
    }

    /* ── Stage 2: derivative at (state + dt/2 * k1) ────────────────────── */
    /*
     * Build a temporary body array displaced by half a step in the k1
     * direction, then evaluate accelerations there.
     */

    for (i = 0; i < n; i++){
        tmp[i]    = b0[i];
        tmp[i].x  = b0[i].x  + 0.5 * dt * k1vx[i];
        tmp[i].y  = b0[i].y  + 0.5 * dt * k1vy[i];
        tmp[i].vx = b0[i].vx + 0.5 * dt * k1ax[i];
        tmp[i].vy = b0[i].vy + 0.5 * dt * k1ay[i];
    }

    compute_accel_array(tmp, n);

    for (i = 0; i < n; i++){
        k2vx[i] = tmp[i].vx;
        k2vy[i] = tmp[i].vy;
        k2ax[i] = tmp[i].ax;
        k2ay[i] = tmp[i].ay;
    }

    /* ── Stage 3: derivative at (state + dt/2 * k2) ────────────────────── */

    for (i = 0; i < n; i++) {
        tmp[i]    = b0[i];
        tmp[i].x  = b0[i].x  + 0.5 * dt * k2vx[i];
        tmp[i].y  = b0[i].y  + 0.5 * dt * k2vy[i];
        tmp[i].vx = b0[i].vx + 0.5 * dt * k2ax[i];
        tmp[i].vy = b0[i].vy + 0.5 * dt * k2ay[i];
    }

    compute_accel_array(tmp, n);

    for (i = 0; i < n; i++) {
        k3vx[i] = tmp[i].vx;
        k3vy[i] = tmp[i].vy;
        k3ax[i] = tmp[i].ax;
        k3ay[i] = tmp[i].ay;
    }

    /* ── Stage 4: derivative at (state + dt * k3) ──────────────────────── */
    for (i = 0; i < n; i++) {
        tmp[i]    = b0[i];
        tmp[i].x  = b0[i].x  + dt * k3vx[i];
        tmp[i].y  = b0[i].y  + dt * k3vy[i];
        tmp[i].vx = b0[i].vx + dt * k3ax[i];
        tmp[i].vy = b0[i].vy + dt * k3ay[i];
    }

    compute_accel_array(tmp, n);

    for (i = 0; i < n; i++) {
        k4vx[i] = tmp[i].vx;
        k4vy[i] = tmp[i].vy;
        k4ax[i] = tmp[i].ax;
        k4ay[i] = tmp[i].ay;
    }

    /* ── Final weighted update ──────────────────────────────────────────── */
    /*
     * RK4 weighted combination:
     *   y(t+dt) = y(t) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
     *
     * We apply this to both positions and velocities.
     * After this update we recompute accelerations so that sim->bodies
     * are in a fully consistent state (positions, velocities, AND
     * accelerations all at time t+dt).
     */

    double sixth_dt = dt / 6.0;
    for (i = 0; i < n; i++){
        sim->bodies[i].x  = b0[i].x
            + sixth_dt * (k1vx[i] + 2.0*k2vx[i] + 2.0*k3vx[i] + k4vx[i]);
        sim->bodies[i].y  = b0[i].y
            + sixth_dt * (k1vy[i] + 2.0*k2vy[i] + 2.0*k3vy[i] + k4vy[i]);
        sim->bodies[i].vx = b0[i].vx
            + sixth_dt * (k1ax[i] + 2.0*k2ax[i] + 2.0*k3ax[i] + k4ax[i]);
        sim->bodies[i].vy = b0[i].vy
            + sixth_dt * (k1ay[i] + 2.0*k2ay[i] + 2.0*k3ay[i] + k4ay[i]);
    }

    /* Recompute accelerations at the final positions so the struct is
     * fully up to date.  This call also serves as the "k1" of the next
     * step if the caller immediately chains another RK4 step.              */

    compute_accel_array(sim->bodies, n);

    /* Advance the simulation clock. */

    sim->t += dt;

}

/* ─────────────────────────────────────────────────────────────────────────────
   Diagnostics
   ─────────────────────────────────────────────────────────────────────────────*/

double nbody_total_energy(const Simulation *sim)
{

    int i, j;
    double energy = 0.0;

    /* ── Kinetic energy: KE = Sum  1/2 m v² ───────────────────────────────── */

    for (i = 0; i < sim->n; i++){
        double v2 = sim->bodies[i].vx * sim->bodies[i].vx
                  + sim->bodies[i].vy * sim->bodies[i].vy;
        energy += 0.5 * sim->bodies[i].mass*v2;
    }

    /* ── Gravitational potential energy: PE = -Σ_{i<j} G*mi*mj / r ─────── */
    /*
     * We use the *un-softened* distance here because potential energy is a
     * diagnostic quantity, not a force computation.  Using the softened
     * distance would make the reported energy inconsistent with what you
     * would calculate analytically from the positions.
     *
     * The potential is negative (bound systems have E < 0 overall).
     */

    for (i = 0; i < sim->n; i++){
        for (j = i +1; j < sim->n; j++){

            double dx = sim->bodies[j].x - sim->bodies[i].x;
            double dy = sim->bodies[j].y - sim->bodies[i].y;
            double dist = sqrt(dx*dx + dy*dy);
            energy -= NBODY_G * sim->bodies[i].mass * sim->bodies[j].mass / dist;

        }
    }

    return energy;
}

void nbody_total_momentum(const Simulation *sim, double *px, double *py)
{
    int i;
    *px = 0.0;
    *py = 0.0;
 
    for (i = 0; i < sim->n; i++) {
        *px += sim->bodies[i].mass * sim->bodies[i].vx;
        *py += sim->bodies[i].mass * sim->bodies[i].vy;
    }
}

void nbody_barycentre(const Simulation *sim, double *cx, double *cy)
{

    int i;
    double total_mass = 0.0;
    *cx = 0.0;
    *cy = 0.0;

    for (i = 0; i < sim->n; i++) {
        double m = sim->bodies[i].mass;
        total_mass += m;
        *cx        += m * sim->bodies[i].x;
        *cy        += m * sim->bodies[i].y;
    }


    /*
     * Guard against a zero-mass simulation (all bodies have mass 0).
     * This should never happen in practice — nbody_add_body does not
     * validate mass — but we avoid a NaN return just in case.
     */
    
    if (total_mass > 0.0){
        *cx /= total_mass;
        *cy /= total_mass;
    }

}