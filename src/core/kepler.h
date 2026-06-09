/*
Code reutilized from https://github.com/dM4XW3LL/Keplerian-Orbital-Simulator
*/

#ifndef KEPLER_H
#define KEPLER_H

/*
This header file contains several functions that might be useful not only in this specific exercise.
Namely the solutions to kepler's equation using both fixed point and newton iterations, as well as conversion from eccentric to true anomaly and computation of orbit position.
*/


/* ── Symbol export macro ────────────────────────────────────────────────────
   Forces every tagged function to be visible in the shared library's dynamic
   symbol table, regardless of compiler visibility flags.
   This is what makes ctypes able to find the functions by name.            */
#if defined(_WIN32) || defined(__CYGWIN__)
    #define KEPLER_API __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
    #define KEPLER_API __attribute__((visibility("default")))
#else
    #define KEPLER_API
#endif
 
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ─────────────────────────────────────────────────────────────────────────────
   Kepler's Equation solvers
   ───────────────────────────────────────────────────────────────────────────*/



/**
 * @brief Solve Kepler's Equation M = E -e*sin(E) via fixed-point iteration.
 * 
 * This function solves Kepler's Equation iteratively using the fixed-point iteration method. In general it needs more iterations than Newton's method to achieve the same accuracy.
 * The equation that is solved iteratively is actually E(i+1) = M + e*sin(E(i)).
 * Returns the eccentric anomaly correspondent to the Mean Anomaly in the input.
 * 
 * @param M Mean Anomaly of the Orbit for a given instant
 * @param e Eccentricity of the Orbit
 * @param iter_count Number of iterations
 * @return double Value of Eccentric Anomaly
 */
double kepler_fixed_point(double M, double e, int *iter_count);


/**
 * @brief Solve Kepler's Equation M = E -e*sin(E) via Newton-Raphson iteration.
 * 
 * This function solves Kepler's Equation iteratively using the Newton-Raphson method. In general it needs less iterations than the fixed point method to achieve the same accuracy.
 * The equation that is solved iteratively is actually E(i+1) = E(i) - g(E(i))/g'(E(i)).
 * Where g(E) = E - e*sin(E)-M and g'(E) = 1 - e*cos(E)
 * Returns the eccentric anomaly correspondent to the Mean Anomaly in the input.
 * 
 * @param M Mean Anomaly of the Orbit for a given instant
 * @param e Eccentricity of the Orbit
 * @param iter_count Number of iterations
 * @return double Value of Eccentric AnomalyCompute distance between two points (x1,y1) and (x2,y2)
 */
double kepler_newton(double M, double e, int *iter_count);


/* ─────────────────────────────────────────────────────────────────────────────
   Anomaly conversions
   ───────────────────────────────────────────────────────────────────────────*/

/**
 * @brief Converts eccentric anomaly E to true anomaly f.
 * 
 * This function converts the calculated value of the Eccentric Anomaly to its correspondent true anomaly.
 * To do this it implements the formula tan(f/2) = sqrt((1+e)/(1-e))tan(E/2)
 * 
 * @param e Eccentricity of the Orbit
 * @param E Eccentric Anomaly
 * @return double True Anomaly
 */
double eccentric_to_true(double e, double E);

/**
 * @brief Combines mean_anomaly_at_time, kepler_newton, and eccentric_to_true to return the true anomaly at time t.
 * 
 * @param t_years Time in years after initial instant
 * @param period  Period of the orbit in years
 * @param M0    Initial Mean Anomaly
 * @param e     Eccentricity of the orbit
 * @param iter_count Iteration count -> Remnant of the code
 * @return double True anomaly at time t
 */
double true_anomaly_at_time(double t_years, double period, double M0, double e, int *iter_count);


/* ─────────────────────────────────────────────────────────────────────────────
   Orbital geometry
   ───────────────────────────────────────────────────────────────────────────*/
/**
 * @brief Compute the orbital radius (distance from focus to body)
 * 
 * Conic section equation: r=a(1-e^2)/(1+e*cos(f))
 * 
 * @param f True anomaly (radians)
 * @param a Semi-major axis (AU)
 * @param e Orbital eccentricity
 * @return double Distance from the focus (star) to the body in AU
 */
double orbit_radius(double f, double a, double e);

/**
 * @brief Computes (x,y) position from true anomaly f, semi-major axis a and eccentricity e.
 * 
 * This function computes the coordinates of an object in a certain orbit knowing its true anomaly and other orbital parameters.
 * The values of x and y are only updated via pointers and should be saved to a file. (There is no return)
 * 
 * @param f True Anomaly
 * @param a Semi-Major Axis
 * @param e Eccentricity
 * @param x Position x = r*cos(f)
 * @param y Position y = r*cos(f)
 */
void orbit_position(double f, double a, double e, double *x, double *y);

/* ─────────────────────────────────────────────────────────────────────────────
   Mean anomaly helpers
   ───────────────────────────────────────────────────────────────────────────*/



/**
 * @brief Computes the mean anomaly for an initial moment based on the orbital parameters of a reference time.
 * 
 * This function determines the mean anomaly at an initial timeframe, knowing the mean anomaly at a certain instant in time afterwards.
 * This initial mean anomaly and its corresponding time will be the "origin" of the referential.
 * 
 * @param lambda Mean Longitude
 * @param phi0 
 * @param years_before_epoch Time in years *AFTER* the reference
 * @param period Period of orbit in years
 * @return double Mean Anomaly at a time t after the reference
 */
double initial_mean_anomaly(double lambda, double phi0, double years_before_epoch, double period);

/**
 * @brief Compute mean anomaly at time t (in years) given initial M0 and period.
 * 
 * This function computes the mean anomaly for any time t after an initial reference in time and its corresponding initial mean anomaly.
 * The default equation for the mean anomaly after some interval t is: M(t) = 2*pi*t/P+M0
 * 
 * @param t_years Time in years after initial instant
 * @param period Period of orbit in years
 * @brief 
 * 
 *
 * @param M0 Initial Mean Anomaly
 * @return double Mean Anomaly after time t_years
 */
double mean_anomaly_at_time(double t_years, double period, double M0);


/* ─────────────────────────────────────────────────────────────────────────────
   Orbital velocity — vis-viva equation
   ───────────────────────────────────────────────────────────────────────────*/

/**
 * @brief Orbital speed in AU/yr via the vis-viva equation.
 * 
 * v=sqrt(GM*(2/r - 1/a))
 * In solar units: GM = 4pi^2*M_star
 * 
 * @param r Current distance from star (AU)
 * @param a Semi-major axis (AU)
 * @param M_Star Stellar mass (Solar Masses)
 * @return double Orbital Speed (in AU/yr)
 */
double orbital_velocity_au_yr(double r, double a,double M_Star);


/**
 * @brief Orbital speed in km/s via the vis-viva equation
 * 
 * Same as orbital_velocity_au_yr, multiplied by the conversion factor
 * 1 AU/yr = 4.74047 km/s.4
 * 
 * @param r Current distance from the star (AU)
 * @param a Semi-major axis (AU)
 * @param M_Star Stellar mass (solar masses)
 * @return double  Orbital speed (km/s)
 */
double orbital_velocity_km_s(double r, double a, double M_Star);

/* ─────────────────────────────────────────────────────────────────────────────
   Orbit progress tracking
   ───────────────────────────────────────────────────────────────────────────*/

/**
 * @brief Time elapsed within the current orbit (years).
 * 
 * Accounts for the initial mean anomaly offset M0, then returns
 * (t - t_offset) mod period, always in [0, period).
 * 
 * @param t_years   Current simulation time (years)
 * @param period    Orbital period (years)
 * @param M0        Initial mean anomaly (radians)
 * @return double   Time into the current orbit (years), in [0, period)
 */
double time_in_current_orbit(double t_years, double period, double M0);


/**
 * @brief Percentage of the current orbit completed.
 * 
 * Returns a value in [0, 100).
 * 0 = at periapsis start, 50 = halfway, approaching 100 = nearly full orbit.
 * 
 * @param t_years    Current simulation time (years)
 * @param period     Orbital period (years)
 * @param M0         Initial mean anomaly (radians)
 * @return double    Orbit completion percentage [0, 100)
 */
double orbit_progress(double t_years, double period, double M0);

/* ─────────────────────────────────────────────────────────────────────────────
   Composite helpers
   ───────────────────────────────────────────────────────────────────────────*/


/**
 * @brief Compute (x, y) position of a planet at time t_years.
 * 
 * This function computes the coordinates (x,y) of a planet in a orbit around the Sun.
 * To achieve this the function utilizes the newton-raphson method internally; iter_count is accumlated.
 * 
 * 
 * @param t_years Time in years after the initial mean anomaly
 * @param a Semi-Major Axis of the Orbit
 * @param e Eccentricity of the Orbit
 * @param period Period of the Orbit in years
 * @param M0 Initial Mean Anomaly
 * @param x x = r*cos(f)
 * @param y y = r*sin(f)
 * @param iter_count This is just necessary for internal reasons, namely because it reuses other functions
 */
void planet_position(double t_years, double a, double e, double period, double M0, double *x, double *y, int *iter_count);



/**
 * @brief Compute distance between two points (x1,y1) and (x2,y2) in two orbits.
 * 
 * This function is able to compute the distance between two points in two different orbits, provided their positions in their respective orbits.
 * However, both orbits should have as the focus/center the same object, for example, the Sun.
 * 
 * @param x1 Coordinate x of first object
 * @param y1 Coordinate y of first object
 * @param x2 Coordinate x of second object
 * @param y2 Coordinate y of second object
 * @return double 
 */
double planet_distance(double x1, double y1, double x2, double y2);


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif




#endif /* KEPLER_H */