/**
 * @file   nbody.h
 * @brief  N-body gravitational simulator — public API.
 *
 * @details
 * This header defines the data structures and functions that make up the
 * gravitational N-body engine.  The engine is self-contained in @c nbody.c
 * and is designed to be compiled as a shared library and called from Python
 * (or any other host) via ctypes / FFI.
 *
 * ### Physical model
 *
 * Every body is treated as a point mass.  The gravitational acceleration
 * acting on body @f$ i @f$ due to all other bodies is:
 *
 * @f[
 *   \mathbf{a}_i = G \sum_{j \neq i} m_j
 *       \frac{\mathbf{r}_j - \mathbf{r}_i}
 *            {\left(|\mathbf{r}_j - \mathbf{r}_i|^2 + \varepsilon^2\right)^{3/2}}
 * @f]
 *
 * where @f$ \varepsilon @f$ (NBODY_SOFTENING) is a small softening length
 * that prevents the denominator from reaching zero during close approaches.
 *
 * ### Unit system
 *
 * All quantities use the same **solar unit system** as kepler.c so that the
 * two libraries interoperate without any conversion:
 *
 * | Quantity  | Unit          | Symbol |
 * |-----------|---------------|--------|
 * | Length    | Astronomical Unit | AU |
 * | Time      | Julian year   | yr     |
 * | Mass      | Solar mass    | M☉     |
 *
 * In these units Newton's gravitational constant is exactly
 * @f$ G = 4\pi^2 \;\mathrm{AU}^3\,\mathrm{yr}^{-2}\,M_\odot^{-1} @f$.
 *
 * ### Integration methods
 *
 * Two integrators are provided:
 *
 * - **Leapfrog (Störmer–Verlet)** — the default for orbital mechanics.
 *   It is a *symplectic* integrator, meaning it conserves a modified
 *   Hamiltonian exactly.  Over long simulations this prevents the secular
 *   (ever-growing) energy drift that plagues non-symplectic methods.
 *
 * - **Runge–Kutta 4 (RK4)** — a classical explicit integrator.  It is
 *   fourth-order accurate (error ∝ dt⁴ per step) and more accurate than
 *   leapfrog for a single large step, but it is *not* symplectic, so energy
 *   drifts slowly over very long runs.  Useful for short, high-accuracy
 *   integrations or validation.
 *
 * ### Typical usage (from C)
 *
 * @code{.c}
 * // 1.  Build a simulation with two bodies.
 * Simulation sim;
 * nbody_init(&sim);
 *
 * // 2.  Add the Sun at rest at the origin.
 * Body sun = {0};
 * sun.mass = 1.0;
 * snprintf(sun.name, NBODY_NAME_LEN, "Sun");
 * nbody_add_body(&sim, &sun);
 *
 * // 3.  Add Earth from orbital elements.
 * Body earth = {0};
 * earth.mass = 3.003e-6;
 * snprintf(earth.name, NBODY_NAME_LEN, "Earth");
 * elements_to_state(1.0, 0.0167, 0.0, 1.0, &earth.x, &earth.y,
 *                   &earth.vx, &earth.vy);
 * nbody_add_body(&sim, &earth);
 *
 * // 4.  Integrate for one year in 500 steps.
 * double dt = 1.0 / 500.0;
 * for (int i = 0; i < 500; i++)
 *     nbody_step_leapfrog(&sim, dt);
 *
 * // 5.  Read back Earth's position.
 * printf("Earth: x=%.6f  y=%.6f AU\n",
 *        sim.bodies[1].x, sim.bodies[1].y);
 *
 * // 6.  Clean up.
 * nbody_free(&sim);
 * @endcode
 *
 * @author  Diogo
 * @date    2026
 * @version 1.0
 */

#ifndef NBODY_H
#define NBODY_H

#include <stddef.h>   /* size_t */

/* ── Symbol export macro ─────────────────────────────────────────────────────
   Ensures every tagged function appears in the shared-library symbol table
   regardless of compiler visibility flags, which is what makes ctypes able
   to look functions up by name at runtime.                                   */
#if defined(_WIN32) || defined(__CYGWIN__)
    #define NBODY_API __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
    #define NBODY_API __attribute__((visibility("default")))
#else
    #define NBODY_API
#endif

/* ── Physical constants ──────────────────────────────────────────────────────*/

/**
 * @brief  Newton's gravitational constant in solar units.
 *
 * @details
 * In the unit system AU / yr / M☉, Kepler's third law gives
 * @f$ P^2 = a^3 / M @f$ directly, which implies
 * @f$ G = 4\pi^2 \approx 39.478 @f$ (AU³ yr⁻² M☉⁻¹).
 *
 * This is the same value used internally by kepler.c so both libraries
 * agree on the same gravitational physics.
 */
#define NBODY_G              39.47841760435743   /* 4 * pi^2 */

/**
 * @brief  Gravitational softening length (AU).
 *
 * @details
 * The softening parameter @f$ \varepsilon @f$ is added in quadrature to
 * the inter-body distance before computing the gravitational force:
 *
 * @f[
 *   r_{\text{soft}} = \sqrt{|\mathbf{r}_j - \mathbf{r}_i|^2 + \varepsilon^2}
 * @f]
 *
 * This prevents the force from diverging when two bodies come extremely
 * close together.  The chosen value (0.001 AU ≈ 150 000 km) is small enough
 * to leave planetary-scale orbits unaffected but large enough to avoid
 * numerical blow-up in typical interactive simulations.
 *
 * For precision research applications this value should be reduced or made
 * user-configurable; for interactive / educational use it is left as a
 * compile-time constant.
 */
#define NBODY_SOFTENING      1e-3                /* AU */

/**
 * @brief  Maximum number of bodies in one simulation.
 *
 * @details
 * Bodies are stored in a fixed-size array inside @ref Simulation to avoid
 * dynamic allocation in the hot integration loop.  64 bodies is far more
 * than any interactive solar-system scenario requires and keeps the struct
 * small enough that copying it is not expensive.
 *
 * Increase this value (and recompile) if you need to simulate denser systems
 * such as star clusters.
 */
#define NBODY_MAX_BODIES     64

/**
 * @brief  Maximum length (including NUL terminator) of a body name string.
 */
#define NBODY_NAME_LEN       64


/* ═══════════════════════════════════════════════════════════════════════════
   Data structures
   ═══════════════════════════════════════════════════════════════════════════*/

/**
 * @brief  Cartesian state vector and physical parameters of one body.
 *
 * @details
 * The N-body integrator works entirely in Cartesian coordinates.  Orbital
 * elements (semi-major axis, eccentricity, …) are a convenient way to
 * *specify* initial conditions, but they are converted to (x, y, vx, vy)
 * before the first integration step via @ref elements_to_state.
 *
 * All physical quantities use the solar unit system defined in the file
 * header.
 *
 * @note  The struct is laid out so that the six doubles (x, y, vx, vy,
 *        ax, ay) are contiguous in memory.  This makes it straightforward
 *        to copy state vectors in bulk if you later want to implement an
 *        RK4 substep that operates on temporary state copies.
 */
typedef struct {
    /* ── Kinematics ─────────────────────────────────────────────────────── */
    double x;       /**< Position x-component (AU).                        */
    double y;       /**< Position y-component (AU).                        */
    double vx;      /**< Velocity x-component (AU yr⁻¹).                  */
    double vy;      /**< Velocity y-component (AU yr⁻¹).                  */
    double ax;      /**< Acceleration x-component (AU yr⁻²).
                         Computed by @ref nbody_compute_accelerations;
                         do not set manually.                               */
    double ay;      /**< Acceleration y-component (AU yr⁻²).
                         Computed by @ref nbody_compute_accelerations;
                         do not set manually.                               */

    /* ── Physical properties ─────────────────────────────────────────────*/
    double mass;    /**< Gravitational mass (M☉).  Must be > 0.           */

    /* ── Metadata ────────────────────────────────────────────────────────*/
    char   name[NBODY_NAME_LEN]; /**< Human-readable label (NUL-terminated). */
} Body;


/**
 * @brief  Container for the entire N-body simulation state.
 *
 * @details
 * A @c Simulation owns all of its bodies and the current simulation time.
 * It is the single object passed to every integrator call, making the API
 * stateless from the caller's perspective — you can serialise and restore
 * a simulation simply by copying this struct.
 *
 * ### Memory ownership
 *
 * Bodies are stored in a fixed-size embedded array (no heap allocation).
 * Call @ref nbody_init before use and @ref nbody_free when done (currently
 * a no-op, but present for forward compatibility if dynamic allocation is
 * added later).
 *
 * ### Thread safety
 *
 * A @c Simulation is **not** thread-safe.  If you need to run multiple
 * simulations concurrently, give each thread its own @c Simulation instance.
 */
typedef struct {
    Body   bodies[NBODY_MAX_BODIES]; /**< Array of all bodies in the system. */
    int    n;                         /**< Number of active bodies (0 … NBODY_MAX_BODIES). */
    double t;                         /**< Current simulation time (yr).  Set to 0 by @ref nbody_init. */
} Simulation;


/* ═══════════════════════════════════════════════════════════════════════════
   Simulation lifecycle
   ═══════════════════════════════════════════════════════════════════════════*/

/**
 * @brief  Initialise a @ref Simulation to a clean, empty state.
 *
 * @details
 * Zero-fills all body slots, sets @c n = 0 and @c t = 0.0.  This must be
 * called before any other function that operates on a @c Simulation.
 *
 * @param[out] sim  Pointer to the Simulation to initialise.  Must not be
 *                  NULL.
 */
NBODY_API void nbody_init(Simulation *sim);


/**
 * @brief  Release any resources held by a @ref Simulation.
 *
 * @details
 * In the current implementation all storage is embedded in the struct
 * (no heap allocations), so this function is effectively a no-op.  It is
 * provided so that callers do not need to change their cleanup code if a
 * future version introduces dynamic allocation (e.g. for >64 bodies).
 *
 * Always call this when you are finished with a simulation.
 *
 * @param[in,out] sim  Pointer to the Simulation to release.
 */
NBODY_API void nbody_free(Simulation *sim);


/**
 * @brief  Add a body to the simulation.
 *
 * @details
 * Copies the body into the next free slot.  The caller does not need to
 * keep the source @c Body alive after this call.
 *
 * The acceleration fields (@c ax, @c ay) in the source body are ignored;
 * they will be computed by the first call to @ref nbody_compute_accelerations.
 *
 * @param[in,out] sim   Simulation to add the body to.
 * @param[in]     body  Body to add (copied by value).
 *
 * @return  0 on success.
 * @return  -1 if the simulation is already full (@c n == NBODY_MAX_BODIES);
 *          the simulation is left unchanged.
 */
NBODY_API int nbody_add_body(Simulation *sim, const Body *body);


/**
 * @brief  Remove a body from the simulation by index.
 *
 * @details
 * The body at position @p idx is removed and the remaining bodies above it
 * are shifted down to fill the gap, preserving their relative order.  The
 * simulation's body count is decremented.
 *
 * @param[in,out] sim  Simulation to modify.
 * @param[in]     idx  Zero-based index of the body to remove
 *                     (must satisfy 0 ≤ idx < sim->n).
 *
 * @return  0 on success.
 * @return  -1 if @p idx is out of range; the simulation is left unchanged.
 */
NBODY_API int nbody_remove_body(Simulation *sim, int idx);


/* ═══════════════════════════════════════════════════════════════════════════
   Initial-condition helpers
   ═══════════════════════════════════════════════════════════════════════════*/

/**
 * @brief  Convert Keplerian orbital elements to a Cartesian state vector.
 *
 * @details
 * This is the bridge between kepler.c and the N-body engine.  When you want
 * to initialise a body on a specific orbit (e.g. "Earth on a 1 AU circular
 * orbit"), you express that orbit in the familiar Keplerian elements and this
 * function converts them to the (x, y, vx, vy) that the integrator needs.
 *
 * ### Derivation
 *
 * Given the eccentric anomaly @f$ E @f$ corresponding to @f$ M_0 @f$
 * (solved via Newton–Raphson, as in kepler.c):
 *
 * @f{align*}{
 *   x   &= a(\cos E - e) \\
 *   y   &= a\sqrt{1-e^2}\sin E
 * @f}
 *
 * Differentiating with respect to time and using @f$ \dot{E} = n/(1-e\cos E) @f$
 * (where @f$ n = 2\pi/P @f$ is the mean motion):
 *
 * @f{align*}{
 *   v_x &= -\frac{na\sin E}{1-e\cos E} \\
 *   v_y &= \frac{na\sqrt{1-e^2}\cos E}{1-e\cos E}
 * @f}
 *
 * This places the body in the orbital plane with the central mass at the
 * *focus* (not the centre) of the ellipse, consistent with kepler.c's
 * coordinate convention.
 *
 * @param[in]  a        Semi-major axis (AU).  Must be > 0.
 * @param[in]  e        Eccentricity.  Must satisfy 0 ≤ e < 1 (elliptic orbit).
 * @param[in]  M0       Initial mean anomaly (radians).  Use 0 for periapsis.
 * @param[in]  M_star   Mass of the central body (M☉).  Used to compute the
 *                      orbital velocity via the vis-viva / Kepler relation;
 *                      must be > 0.
 * @param[out] x        Initial x-position (AU).
 * @param[out] y        Initial y-position (AU).
 * @param[out] vx       Initial x-velocity (AU yr⁻¹).
 * @param[out] vy       Initial y-velocity (AU yr⁻¹).
 */
NBODY_API void elements_to_state(double a, double e, double M0,
                                  double M_star,
                                  double *x,  double *y,
                                  double *vx, double *vy);


/* ═══════════════════════════════════════════════════════════════════════════
   Core physics
   ═══════════════════════════════════════════════════════════════════════════*/

/**
 * @brief  Compute the gravitational acceleration on every body.
 *
 * @details
 * For each pair (i, j) with i ≠ j, the Newtonian gravitational acceleration
 * on body @f$ i @f$ due to body @f$ j @f$ is:
 *
 * @f[
 *   \mathbf{a}_{ij} = G m_j
 *       \frac{\mathbf{r}_j - \mathbf{r}_i}
 *            {\left(|\mathbf{r}_j - \mathbf{r}_i|^2 + \varepsilon^2\right)^{3/2}}
 * @f]
 *
 * The total acceleration on body @f$ i @f$ is the vector sum over all
 * @f$ j \neq i @f$.  Because the force is symmetric (Newton's third law),
 * the inner loop only evaluates each pair once and adds the result to both
 * bodies simultaneously, halving the number of distance computations.
 *
 * The result is written directly into each body's @c ax / @c ay fields.
 *
 * **Complexity:** O(N²) — unavoidable for the exact direct sum.  For
 * N ≤ 64 (the current limit) this is very fast; for larger N a Barnes–Hut
 * tree (O(N log N)) would be needed.
 *
 * @param[in,out] sim  Simulation whose body accelerations are updated.
 */
NBODY_API void nbody_compute_accelerations(Simulation *sim);


/* ═══════════════════════════════════════════════════════════════════════════
   Integrators
   ═══════════════════════════════════════════════════════════════════════════*/

/**
 * @brief  Advance the simulation by one timestep using the Leapfrog integrator.
 *
 * @details
 * The leapfrog (Störmer–Verlet) method is the recommended integrator for
 * orbital mechanics because it is **symplectic**: it conserves a modified
 * Hamiltonian exactly, which means total energy oscillates around a constant
 * rather than drifting upward or downward over long simulations.
 *
 * ### Algorithm (velocity-Verlet form of leapfrog)
 *
 * Each call performs three sub-steps:
 *
 * 1. **Half velocity kick** — advance velocities by half a timestep using
 *    the accelerations already stored in each body:
 *    @f[ v_{i} \leftarrow v_{i} + \tfrac{1}{2} a_{i}(t)\,\Delta t @f]
 *
 * 2. **Full position drift** — advance positions by a full timestep using
 *    the half-kicked velocities:
 *    @f[ r_{i} \leftarrow r_{i} + v_{i}\,\Delta t @f]
 *
 * 3. **Recompute accelerations** — evaluate the gravitational forces at the
 *    new positions via @ref nbody_compute_accelerations.
 *
 * 4. **Second half velocity kick** — advance velocities by another half step
 *    using the *new* accelerations:
 *    @f[ v_{i} \leftarrow v_{i} + \tfrac{1}{2} a_{i}(t+\Delta t)\,\Delta t @f]
 *
 * On the very first call, if @c ax / @c ay are still zero (i.e. they have
 * not been computed yet), the first half-kick is a no-op and the function
 * naturally bootstraps correctly.  If you want to be explicit, call
 * @ref nbody_compute_accelerations once before the first step.
 *
 * ### Choosing a timestep
 *
 * As a rule of thumb, the timestep should be less than ~1% of the shortest
 * orbital period in the system.  For the inner Solar System (Mercury at
 * 0.24 yr) this means @f$ \Delta t \lesssim 0.002 \;\mathrm{yr} @f$.  Using
 * a larger step will cause orbits to precess artificially or become unstable.
 *
 * @param[in,out] sim  Simulation to advance.
 * @param[in]     dt   Timestep (yr).  Must be > 0.
 */
NBODY_API void nbody_step_leapfrog(Simulation *sim, double dt);


/**
 * @brief  Advance the simulation by one timestep using the RK4 integrator.
 *
 * @details
 * The classical 4th-order Runge–Kutta method evaluates the equations of
 * motion at four intermediate points and combines them with weights that
 * cancel error terms up to @f$ O(\Delta t^4) @f$.  It is more accurate than
 * leapfrog for a single large step, but it is **not symplectic**: over very
 * long simulations, energy drifts slowly (typically upward, as if the system
 * is slowly gaining kinetic energy).
 *
 * ### When to use RK4 vs Leapfrog
 *
 * | Scenario                               | Recommended integrator |
 * |----------------------------------------|------------------------|
 * | Long-term orbital stability             | Leapfrog               |
 * | Short high-precision integration        | RK4                    |
 * | Validation / cross-checking leapfrog   | RK4                    |
 * | Highly eccentric orbits (large dt)     | RK4 (more stable)      |
 *
 * ### Algorithm
 *
 * Let the state vector be @f$ \mathbf{y} = (\mathbf{r}, \mathbf{v}) @f$
 * and the derivative @f$ \dot{\mathbf{y}} = (\mathbf{v}, \mathbf{a}(\mathbf{r})) @f$.
 * RK4 computes four derivative evaluations:
 *
 * @f{align*}{
 *   k_1 &= \dot{\mathbf{y}}(\mathbf{y}) \\
 *   k_2 &= \dot{\mathbf{y}}(\mathbf{y} + \tfrac{\Delta t}{2} k_1) \\
 *   k_3 &= \dot{\mathbf{y}}(\mathbf{y} + \tfrac{\Delta t}{2} k_2) \\
 *   k_4 &= \dot{\mathbf{y}}(\mathbf{y} + \Delta t\,k_3)
 * @f}
 *
 * and then advances the state:
 *
 * @f[
 *   \mathbf{y}(t + \Delta t) = \mathbf{y}(t)
 *     + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)
 * @f]
 *
 * In this implementation, each @f$ k @f$ evaluation allocates a temporary
 * copy of the full body array on the stack, so it does not modify @c sim
 * until the final weighted update.
 *
 * @param[in,out] sim  Simulation to advance.
 * @param[in]     dt   Timestep (yr).  Must be > 0.
 */
NBODY_API void nbody_step_rk4(Simulation *sim, double dt);


/* ═══════════════════════════════════════════════════════════════════════════
   Diagnostics
   ═══════════════════════════════════════════════════════════════════════════*/

/**
 * @brief  Compute the total mechanical energy of the system.
 *
 * @details
 * The total energy is the sum of kinetic and gravitational potential energy:
 *
 * @f[
 *   E = \sum_i \tfrac{1}{2} m_i v_i^2
 *     - G \sum_{i < j} \frac{m_i m_j}{|\mathbf{r}_j - \mathbf{r}_i|}
 * @f]
 *
 * For a closed, isolated system integrated with a symplectic method (leapfrog)
 * this quantity should oscillate around a constant.  A secular drift indicates
 * either too large a timestep or numerical precision issues.
 *
 * **Usage tip:** Call this function before and after a long integration run
 * and compare the values.  A relative change of less than 0.01% over a full
 * orbital period is generally considered acceptable.
 *
 * @param[in] sim  Simulation to analyse (read-only).
 * @return         Total mechanical energy (M☉ AU² yr⁻²).
 */
NBODY_API double nbody_total_energy(const Simulation *sim);


/**
 * @brief  Compute the total linear momentum of the system.
 *
 * @details
 * The total momentum is:
 * @f[
 *   \mathbf{p} = \sum_i m_i \mathbf{v}_i
 * @f]
 *
 * For an isolated system this should be conserved to machine precision by
 * any integrator (it is not related to symplecticity).  A drift in momentum
 * indicates a bug in the force computation.
 *
 * @param[in]  sim  Simulation to analyse (read-only).
 * @param[out] px   x-component of total momentum (M☉ AU yr⁻¹).
 * @param[out] py   y-component of total momentum (M☉ AU yr⁻¹).
 */
NBODY_API void nbody_total_momentum(const Simulation *sim,
                                     double *px, double *py);


/**
 * @brief  Compute the kinetic energy of a single body.
 *
 * @details
 * @f[
 *   KE_i = \tfrac{1}{2} m_i \left( v_{x,i}^2 + v_{y,i}^2 \right)
 * @f]
 *
 * This is provided as a convenience for GUIs and analysis tools that want
 * to display per-body energy without having to recompute the full system
 * energy via @ref nbody_total_energy.  It does not include the gravitational
 * potential energy contributed by this body's interaction with others.
 *
 * @param[in] sim  Simulation containing the body (read-only).
 * @param[in] idx  Zero-based index of the body (must satisfy 0 ≤ idx < sim->n).
 *
 * @return  Kinetic energy of body @p idx (M☉ AU² yr⁻²).
 * @return  0.0 if @p idx is out of range.
 */
NBODY_API double nbody_body_kinetic_energy(const Simulation *sim, int idx);


/**
 * @brief  Compute the position of the system's barycentre (centre of mass).
 *
 * @details
 * @f[
 *   \mathbf{r}_{\text{cm}} = \frac{\sum_i m_i \mathbf{r}_i}{\sum_i m_i}
 * @f]
 *
 * In a true solar system the barycentre lies just outside the surface of
 * the Sun due to Jupiter's mass.  Tracking it is useful for:
 * - Correcting coordinate frames (display bodies relative to barycentre).
 * - Verifying that the integrator conserves the barycentre position when
 *   the total momentum is non-zero.
 *
 * @param[in]  sim  Simulation to analyse (read-only).
 * @param[out] cx   x-coordinate of the barycentre (AU).
 * @param[out] cy   y-coordinate of the barycentre (AU).
 */
NBODY_API void nbody_barycentre(const Simulation *sim,
                                 double *cx, double *cy);


#endif /* NBODY_H */