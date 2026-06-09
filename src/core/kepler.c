/*
Code reutilized from https://github.com/dM4XW3LL/Keplerian-Orbital-Simulator
*/

#include <math.h>
#include "kepler.h"

#define TOL 1e-9


/*
We define the tolerance of the program so that the accuracy of the iteration methods is smaller than this number.
Since we need AT LEAST 9 decimal places, we will need to use doubles.   
*/

/* ─────────────────────────────────────────────────────────────────────────────
   Kepler's Equation solvers
   ───────────────────────────────────────────────────────────────────────────*/


double kepler_fixed_point(double M, double e, int *iter_count){

    /*Initial Guess: E_0 = M for e<0.8 else pi*/
    double E; //Declares variable of initial guess of eccentric anomaly
    if(e<0.8){
        E = M; 
    }
    else{
        E = M_PI;
    }

    double E_new; //Declares variable for next iteration step

    while(1){

        E_new = M+e*sin(E); //Fixed Point Iteration Algorithm
        (*iter_count)++; //Iteration Counter
        if (fabs(E_new-E)<=TOL){ //Condition for break
            break;
        }

        E = E_new; //If not within accuracy range, do another step

    }

    return E_new;
}

double kepler_newton(double M, double e, int *iter_counter){

    /*Initial Guess: E_0 = M for e<0.8 else pi*/
    double E;
    if(e<0.8){
        E = M;
    }
    else{
        E = M_PI;
    }

    double E_new, g, g_prime; //Declare variables for next step value and g and g' values for E in each step

    while (1){

        g = E - e*sin(E)-M;
        g_prime = 1.0-e*cos(E);
        E_new = E - g/g_prime;
        (*iter_counter)++;
        if(fabs(E_new-E)<=TOL){
            break;
        }
        E = E_new;

    }

    return E_new;

}
/* ─────────────────────────────────────────────────────────────────────────────
   Anomaly conversions
   ───────────────────────────────────────────────────────────────────────────*/


double eccentric_to_true(double e, double E){

    /* THe standard formula to convert from Eccentric to true anomaly is:
    tan(f/2) = sqrt((1+e)/(1-e))*tan(E/2) 

    However, this equation has a problem... the ambiguity of the quadrant since atan(x) only returns values from -pi/2 to pi/2, so we would only get half orbits.
    This can be solved splitting the tan(E/2) into sin and cos components and using the atan2(X,Y) function which already takes quadrants into consideration

    
    
    */
    
    double sin_term = sqrt(1.0+e)*sin(E/2.0);
    double cos_term = sqrt(1.0-e)*cos(E/2.0);
    return 2.0*atan2(sin_term,cos_term);

}

double true_anomaly_at_time(double t_years, double period, double M0, double e, int *iter_count){

        /*
     * Full pipeline: t → M → E → f
     * Combines mean_anomaly_at_time, kepler_newton, and eccentric_to_true
     * into one call for convenience.
     * Returns the true anomaly in radians at time t_years.
     */
    double M = mean_anomaly_at_time(t_years,period,M0);
    double E = kepler_newton(M,e,iter_count);
    return eccentric_to_true(e,E);

}

/* ─────────────────────────────────────────────────────────────────────────────
   Orbital geometry
   ───────────────────────────────────────────────────────────────────────────*/

double orbit_radius(double f, double a, double e){

    /*
     * Distance from the focus (star) to the orbiting body.
     * Conic section equation: r = a(1 - e²) / (1 + e·cos(f))
     */
    return a*(1.0-e*e)/(1.0+e*cos(f));

}

void orbit_position(double f, double a, double e, double *x, double *y){


    double r = a*(1.0-e*e)/(1.0+e*cos(f));
    *x = r*cos(f);
    *y = r*sin(f);
}

/* ─────────────────────────────────────────────────────────────────────────────
   Mean anomaly helpers
   ───────────────────────────────────────────────────────────────────────────*/


double initial_mean_anomaly(double lambda, double phi0, double years_before_epoch, double period){

    return lambda-phi0-years_before_epoch*2.0*M_PI/period;
}

double mean_anomaly_at_time(double t_years, double period, double M0){

    return 2.0*M_PI*t_years/period+M0;

}

/* ─────────────────────────────────────────────────────────────────────────────
   Orbital velocity — vis-viva equation
   ───────────────────────────────────────────────────────────────────────────*/

double orbital_velocity_au_yr(double r, double a,double M_Star){

    /*
     * Vis-viva equation: v² = GM · (2/r - 1/a)
     *
     * In solar units (AU, years, solar masses):
     *   GM = 4π² · M_star   [AU³/yr²]
     *
     * Returns speed in AU/yr.
     */
    double GM = 4.0*M_PI*M_PI*M_Star;
    return sqrt(GM*(2.0/r-1.0/a));

}

double orbital_velocity_km_s(double r, double a, double M_Star){

    return orbital_velocity_au_yr(r,a,M_Star)*4.74047;
}

/* ─────────────────────────────────────────────────────────────────────────────
   Orbit progress tracking
   ───────────────────────────────────────────────────────────────────────────*/

double time_in_current_orbit(double t_years, double period, double M0)
{
    /*
     * How far into its current orbit is the body, in years?
     * We subtract the time offset implied by M0 (M0 / 2π · period),
     * then take the remainder modulo the period.
     * Result is always in [0, period).
     */
    double t_offset = M0 / (2.0 * M_PI) * period;
    double t_adj    = t_years - t_offset;
 
    /* fmod can return negative values if t_adj < 0, so we correct for that */
    double t_mod = fmod(t_adj, period);
    if (t_mod < 0.0)
        t_mod += period;
 
    return t_mod;
}

double orbit_progress(double t_years, double period, double M0)
{
    /*
     * Percentage of the current orbit completed, in [0, 100).
     */
    return 100.0 * time_in_current_orbit(t_years, period, M0) / period;
}

/* ─────────────────────────────────────────────────────────────────────────────
   Composite helpers
   ───────────────────────────────────────────────────────────────────────────*/


void planet_position(double t_years, double a, double e, double period, double M0, double *x, double *y, int *iter_count){

    double M = mean_anomaly_at_time(t_years,period,M0);
    double E = kepler_newton(M, e, iter_count);
    double f = eccentric_to_true(e,E);
    orbit_position(f,a,e,x,y);

}

double planet_distance(double x1,double y1,double x2,double y2){

    double dx = x2-x1;
    double dy = y2-y1;
    return sqrt(dx*dx+dy*dy);

}