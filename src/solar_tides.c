/** * @file central_force.c
 * @brief   A general central force.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Central Force$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                None
 * C Example               :ref:`c_example_central_force`
 * Python Example          `CentralForce.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CentralForce.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds a general central acceleration of the form a=Acentral*r^gammacentral, outward along the direction from a central particle to the body.
 * Effect is turned on by adding Acentral and gammacentral parameters to a particle, which will act as the central body for the effect,
 * and will act on all other particles.
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * Acentral (double)             Yes         Normalization for central acceleration.
 * gammacentral (double)         Yes         Power index for central acceleration.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_solar_tides(struct reb_simulation* const sim, struct reb_particle* const particles){
    const int N = sim->N;
    const int source_index = N-1;
    const struct reb_particle source = particles[source_index];
    //calculate center of mass excluding sun particle
    struct reb_particle com = reb_simulation_com_range(sim,0,N-1); 
    // struct reb_particle com = reb_get_com(sim);
    // calculate acceleration felt by COM
    const double dx = com.x - source.x;
    const double dy = com.y - source.y;
    const double dz = com.z - source.z;
    const double r = sqrt(dx*dx + dy*dy + dz*dz);
    const double prefac = sim->G*source.m/pow(r,3);
    const double ax_com = prefac*dx;
    const double ay_com = prefac*dy;
    const double az_com = prefac*dz;

    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r = sqrt(dx*dx + dy*dy + dz*dz);
        const double prefac = sim->G*source.m/pow(r,3);

        particles[i].ax -= (prefac*dx - ax_com);
        particles[i].ay -= (prefac*dy - ay_com);
        particles[i].az -= (prefac*dz - az_com);
        // particles[source_index].ax -= p.m/source.m*prefac*dx; #doesnt matter b/c we are deleting
        // particles[source_index].ay -= p.m/source.m*prefac*dy;
        // particles[source_index].az -= p.m/source.m*prefac*dz;
    }
}

void rebx_solar_tides(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles){
    /* default values */
    double* m_sun = rebx_get_param(sim->extras, force->ap, "m_sun");
    if (m_sun == NULL){
        reb_simulation_error(sim, "Need to specify the a mass\n");
    }
    double* a_sun= rebx_get_param(sim->extras, force->ap, "a_sun");
    if (a_sun == NULL){
        reb_simulation_error(sim, "Need to specify the a semimajor axis\n");
    }
    double* e_sun= rebx_get_param(sim->extras, force->ap, "e_sun");
    if (e_sun== NULL){
        reb_simulation_error(sim, "Need to specify an eccentricity\n");
    }
    double* i_sun= rebx_get_param(sim->extras, force->ap, "i_sun");
    if (i_sun == NULL){
        reb_simulation_error(sim, "Need to specify an inclination.\n");
    }
    double* Omega_sun= rebx_get_param(sim->extras, force->ap, "Omega_sun");
    if (Omega_sun== NULL){
        reb_simulation_error(sim, "Need to specify a Longitude of ascending node\n");
    }
    double* omega_sun= rebx_get_param(sim->extras, force->ap, "omega_sun");
    if (omega_sun == NULL){
        reb_simulation_error(sim, "Need specify an argument of periapse.\n");
    }
    double* T_peri= rebx_get_param(sim->extras, force->ap, "T_peri");
    if (T_peri == NULL){
        reb_simulation_error(sim, "Need to specify a time of pericenter passage.\n");
    }
    // printf("%lf\n", *m_sun);
    // printf("%lf\n", a_sun);
    // printf("%lf\n", e_sun);
    // printf("%lf\n", i_sun);

    double n = sqrt(sim->G * (*m_sun) / pow(*a_sun,3));
    // printf("n=%lf t=%lf T=%lf\n", 1e6*n, sim->t, *T_peri);

    double M = n * (sim->t - (*T_peri));
    // printf("M=%lf M=%lf\n", M, *T_peri);

    /* 
    hack: add sun particle simulation. this allows us to easily get x,y,z position relative to system center
    then we remove the sun from the simulation.
    */
    // add our sun
    // reb_simulation_add_fmt(sim, "m a e inc Omega omega T", *m_sun, *a_sun, *e_sun, *i_sun, *Omega_sun, *omega_sun, *T_peri);
    reb_simulation_add_fmt(sim, "m a e inc Omega omega M", *m_sun, *a_sun, *e_sun, *i_sun, *Omega_sun, *omega_sun, M);
    
    // struct reb_orbit o =  reb_orbit_from_particle(sim->G, sim->particles[sim->N-1], reb_simulation_com_range(sim,0,sim->N-1));
    // printf("M=%f f=%f\n", o.M, n);
    // struct reb_particle reb_particle_from_fmt
    //compute accelerations of particles felt due to sun RELATIVE to COM 
    rebx_calculate_solar_tides(sim, particles); // only calculates force if a particle has both Acentral and gammacentral parameters set.

    //remove
    reb_simulation_remove_particle(sim, sim->N-1, 1);
    // for (int i=0; i<N; i++){
    //     const double* const Acentral = rebx_get_param(sim->extras, particles[i].ap, "Acentral");
    //     if (Acentral != NULL){
    //         const double* const gammacentral = rebx_get_param(sim->extras, particles[i].ap, "gammacentral");
    //         if (gammacentral != NULL){
    //             rebx_calculate_central_force(sim, particles); // only calculates force if a particle has both Acentral and gammacentral parameters set.
    //         }
    //     }
    // }
}

// static double rebx_calculate_central_force_potential(struct reb_simulation* const sim, const double A, const double gamma, const int source_index){
//     const struct reb_particle* const particles = sim->particles;
// 	const int _N_real = sim->N - sim->N_var;
//     const struct reb_particle source = particles[source_index];
//     double H = 0.;
// 	for (int i=0;i<_N_real;i++){
// 		if(i == source_index){
//             continue;
//         }
//         const struct reb_particle p = particles[i];
//         const double dx = p.x - source.x;
//         const double dy = p.y - source.y;
//         const double dz = p.z - source.z;
//         const double r2 = dx*dx + dy*dy + dz*dz;

//         if (fabs(gamma+1.) < DBL_EPSILON){ // F propto 1/r
//             H -= p.m*A*log(sqrt(r2));
//         }
//         else{
//             H -= p.m*A*pow(r2, (gamma+1.)/2.)/(gamma+1.);
//         }
//     }		
//     return H;
// }

// double rebx_central_force_potential(struct rebx_extras* const rebx){
//     if (rebx->sim == NULL){
//         rebx_error(rebx, ""); // rebx_error gives meaningful err
//         return 0;
//     }
//     struct reb_simulation* sim = rebx->sim;
//     const int N_real = sim->N - sim->N_var;
//     struct reb_particle* const particles = sim->particles;
//     double Htot = 0.;
//     for (int i=0; i<N_real; i++){
//         const double* const Acentral = rebx_get_param(rebx, particles[i].ap, "Acentral");
//         if (Acentral != NULL){
//             const double* const gammacentral = rebx_get_param(rebx, particles[i].ap, "gammacentral");
//             if (gammacentral != NULL){
//                 Htot += rebx_calculate_central_force_potential(sim, *Acentral, *gammacentral, i);
//             }
//         }
//     }
//     return Htot;
// }

// double rebx_central_force_Acentral(const struct reb_particle p, const struct reb_particle primary, const double pomegadot, const double gamma){
//     struct reb_simulation* sim = p.sim;
//     const double G = sim->G;
//     const struct reb_orbit o = reb_tools_particle_to_orbit(G, p, primary);
//     if (fabs(gamma+2.) < DBL_EPSILON){  // precession goes to 0 at r^-2, so A diverges for gamma=-2
//         reb_simulation_error(sim, "Precession vanishes for force law varying as r^-2, so can't initialize Acentral from a precession rate for gamma=-2)\n");
//         return 0.;
//     }
//     return G*primary.m*pomegadot/(1.+gamma/2.)/pow(o.d, gamma+2.)/o.n;
// }
