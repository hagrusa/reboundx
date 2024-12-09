/** * @file gravity_second_order.c
 * @brief   Add C20 and C22 gravity harmonics to particles
 * @author  Harrison Agrusa <hagrusa@oca.eu>
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
 * $Gravity Fields$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
 * Based on                None
 * C Example               :ref:`c_example_C20`
 * Python Example          `C20.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/C20.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds 2nd degree & order gravity terms (C20 and C22) to particle. Current implementation for C22 requires a rotating frame such that primary body is always in principle axis alignment
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
 * C20 (double)                  No          C20 coefficient (equivalent to -C20)
 * C22 (double)                  No          C22 coefficient 
 * R_eq (double)                 No          Equatorial radius of nonspherical body used for normalizing harmonics
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_C20_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double C20, const double R_eq, const int source_index){
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        //The following block is the original code from Dan Tamayo's J2 example
        // const struct reb_particle p = particles[i];
        // const double dx = p.x - source.x;
        // const double dy = p.y - source.y;
        // const double dz = p.z - source.z;
        // const double r2 = dx*dx + dy*dy + dz*dz;
        // const double r = sqrt(r2);
        // const double costheta2 = dz*dz/r2;
        // const double prefac = 3.*C20*R_eq*R_eq/r2/r2/r/2.;
        // const double fac = 5.*costheta2-1.;

        // particles[i].ax += G*source.m*prefac*fac*dx;
        // particles[i].ay += G*source.m*prefac*fac*dy;
        // particles[i].az += G*source.m*prefac*(fac-2.)*dz;
        // particles[source_index].ax -= G*p.m*prefac*fac*dx;
        // particles[source_index].ay -= G*p.m*prefac*fac*dy;
        // particles[source_index].az -= G*p.m*prefac*(fac-2.)*dz;

        const struct reb_particle p = particles[i];
        const double x = p.x - source.x;
        const double y = p.y - source.y;
        const double z = p.z - source.z;
        const double r2 = x*x + y*y + z*z;
        const double r = sqrt(r2);
        const double prefac1 = G*source.m*C20*R_eq*R_eq/r2/r2/r;
        const double prefac2 = 5*(r2 - 3*z*z)/2/r2;
        particles[i].ax += prefac1*x*(prefac2-1);
        particles[i].ay += prefac1*y*(prefac2-1);
        particles[i].az += prefac1*z*(prefac2+2);
        particles[source_index].ax -= (p.m/source.m)*prefac1*x*(prefac2-1);
        particles[source_index].ay -= (p.m/source.m)*prefac1*y*(prefac2-1);
        particles[source_index].az -= (p.m/source.m)*prefac1*z*(prefac2+2);
        }
}

static void rebx_C20(struct rebx_extras* const rebx, struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const C20 = rebx_get_param(rebx, particles[i].ap, "C20");
        if (C20 != NULL){
            const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
            if (R_eq != NULL){
                rebx_calculate_C20_force(sim, particles, N, *C20, *R_eq,i); 
            }
        }
    }
}

static void rebx_calculate_C22_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double C22, const double R_eq, const double wPrim, const int source_index){
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        // const struct reb_particle p = particles[i];
        // const double dx = p.x - source.x;
        // const double dy = p.y - source.y;
        // const double dz = p.z - source.z;
        // const double r2 = dx*dx + dy*dy + dz*dz;
        // const double r = sqrt(r2);
        // const double costheta2 = dz*dz/r2;
        // const double prefac = 5.*C22*R_eq*R_eq*R_eq*R_eq/r2/r2/r2/r/8.;
        // const double fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;

        // particles[i].ax += G*source.m*prefac*fac*dx;
        // particles[i].ay += G*source.m*prefac*fac*dy;
        // particles[i].az += G*source.m*prefac*(fac+12.-28.*costheta2)*dz;
        // particles[source_index].ax -= G*p.m*prefac*fac*dx;
        // particles[source_index].ay -= G*p.m*prefac*fac*dy;
        // particles[source_index].az -= G*p.m*prefac*(fac+12.-28.*costheta2)*dz


        // lets calculate everything in the primary's body fixed frame
        // then rotate everything to the inertial frame
        const double t = sim->t;
        const double phi = wPrim * t;
        const double sinphi = sin(phi);
        const double cosphi = cos(phi);

        const struct reb_particle p = particles[i];
        const double x_i = p.x - source.x;
        const double y_i = p.y - source.y;
        const double z_i = p.z - source.z;
        //rotate to body-fixed frame:
        const double x = x_i*cosphi + y_i*sinphi; //body fixed frame
        const double y = -x_i*sinphi + y_i*cosphi;
        const double z = z_i;
        
        const double r2 = x*x + y*y + z*z;
        const double r = sqrt(r2);
        const double prefac1 = G*source.m*C22*R_eq*R_eq/r2/r2/r;
        const double prefac2 = 15*(x*x-y*y)/r2;

        const double ax = -prefac1*x*(prefac2-6.0);
        const double ay = -prefac1*y*(prefac2+6.0);
        const double az = -prefac1*prefac2*z;

        const double ax_i = ax*cosphi - ay*sinphi; //inertial acceleration
        const double ay_i = ax*sinphi + ay*cosphi;
        const double az_i = az;
        

        particles[i].ax += ax_i;
        particles[i].ay += ay_i;
        particles[i].az += az_i;
        particles[source_index].ax -= (p.m/source.m)*ax_i;
        particles[source_index].ay -= (p.m/source.m)*ay_i;
        particles[source_index].az -= (p.m/source.m)*az_i;

    }
}

static void rebx_C22(struct rebx_extras* const rebx, struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const C22 = rebx_get_param(rebx, particles[i].ap, "C22");
        if (C22 != NULL){
            const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
            if (R_eq != NULL){
                const double* const wPrim = rebx_get_param(rebx, particles[i].ap, "wPrim");
                if (wPrim != NULL){
                    rebx_calculate_C22_force(sim, particles, N, *C22, *R_eq, *wPrim, i); 
                }
            }
        }
    }
}

void rebx_gravity_second_order(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    rebx_C20(sim->extras, sim, force, particles, N);
    rebx_C22(sim->extras, sim, force, particles, N);
}

static double rebx_calculate_C20_potential(struct reb_simulation* const sim, const double C20, const double R_eq, const int source_index){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    double H = 0.;
	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        // const struct reb_particle p = particles[i];
        // const double dx = p.x - source.x;
        // const double dy = p.y - source.y;
        // const double dz = p.z - source.z;
        // const double r2 = dx*dx + dy*dy + dz*dz;
        // const double r = sqrt(r2);
        // const double costheta2 = dz*dz/r2;
        // const double prefac = G*p.m*source.m*R_eq*R_eq/r2/r*C20;
        // const double P2 = 0.5*(3.*costheta2-1.);
        // H += prefac*P2;

        const struct reb_particle p = particles[i];
        const double x = p.x - source.x;
        const double y = p.y - source.y;
        const double z = p.z - source.z;
        const double r2 = x*x + y*y + z*z;
        const double r = sqrt(r2);
        // const double costheta2 = dz*dz/r2;
        // const double prefac = G*p.m*source.m*R_eq*R_eq/r2/r*C20;
        // const double P2 = 0.5*(3.*costheta2-1.);
        const double prefac = 0.5*G*source.m*C20*R_eq*R_eq/r2/r2/r;
        H += prefac*(r2 - 3*z*z);
    }		
    return H;
}

static double rebx_C20_potential(struct rebx_extras* const rebx, struct reb_simulation* const sim){
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=0; i<N_real; i++){
        const double* const C20 = rebx_get_param(rebx, particles[i].ap, "C20");
        if (C20 != NULL){
            const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
            if (R_eq != NULL){
                Htot += rebx_calculate_C20_potential(sim, *C20, *R_eq, i);
            }
        }
    }
    return Htot;
}

static double rebx_calculate_C22_potential(struct reb_simulation* const sim, const double C22, const double R_eq, const int source_index){
    // THIS FUNCTION IS WRONG
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    double H = 0.;
	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        // const struct reb_particle p = particles[i];
        // const double dx = p.x - source.x;
        // const double dy = p.y - source.y;
        // const double dz = p.z - source.z;
        // const double r2 = dx*dx + dy*dy + dz*dz;
        // const double r = sqrt(r2);
        // const double costheta2 = dz*dz/r2;
        // const double prefac = G*p.m*source.m*R_eq*R_eq*R_eq*R_eq/r2/r2/r*C22;
        // const double P4 = (35.*costheta2*costheta2 - 30.*costheta2+3.)/8.;
        // H += prefac*P4;
        const struct reb_particle p = particles[i];
        const double x = p.x - source.x;
        const double y = p.y - source.y;
        const double z = p.z - source.z;
        const double r2 = x*x + y*y + z*z;
        const double r = sqrt(r2);
        const double prefac = -3*G*source.m*C22*R_eq*R_eq/r2/r2/r;
        H += prefac*(x*x-y*y);
    }		
    return H;
}

static double rebx_C22_potential(struct rebx_extras* const rebx, struct reb_simulation* const sim){
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=0; i<N_real; i++){
        const double* const C22 = rebx_get_param(rebx, particles[i].ap, "C22");
        if (C22 != NULL){
            const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
            if (R_eq != NULL){
                Htot += rebx_calculate_C22_potential(sim, *C22, *R_eq, i);
            }
        }
    }
    return Htot;
}

double rebx_second_order_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    double H = rebx_C20_potential(rebx, rebx->sim);
    H += rebx_C22_potential(rebx, rebx->sim);
    return H;
}
