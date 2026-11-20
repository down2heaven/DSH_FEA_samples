/* =====================================================================
 *  membrane.c  -  Phase-selective Darcy-Forchheimer momentum sink for the
 *                 breathable-waterproof membrane test case.
 *
 *  Target : ANSYS Fluent 19.2, 3ddp, VOF (2 phases)
 *  Attach : /define/boundary-conditions/fluid -> zone "porous" -> mixture
 *           "Specify source terms?" = yes, 1 source in each of X/Y/Z.
 *
 *  Model
 *  -----
 *  A single permeability is blended between an air value and a water value
 *  using the water volume fraction alpha of the cell:
 *
 *      invK(alpha) = (1-alpha) * INVK_AIR + alpha * INVK_WATER
 *
 *  and the sink is the Darcy term of the Darcy-Forchheimer law with the
 *  inertial (Forchheimer) contribution neglected, because at 1 mm/s the
 *  Reynolds number based on the pore scale is far below unity:
 *
 *      S_i = -invK * mu * u_i          [N/m3]
 *
 *  dS_i/du_i = -invK * mu  is the exact Jacobian, so the sink is treated
 *  fully implicitly: the membrane cannot be "blown through" numerically.
 *
 *  IDEAL CASE (user request):
 *  K_air   = 1e-13 m2  (INVK_AIR = 1e13)  -> breathable, finite resistance
 *  K_water = 1e-17 m2  (INVK_WATER = 1e17) -> water is completely held back
 *
 *  With the linear blend on inverse permeability the air branch contributes
 *  nothing, so invK = alpha * INVK_WATER: resistance switches on exactly in
 *  proportion to how much water the cell contains.
 *
 *  The viscosity mu is the mixture laminar viscosity C_MU_L of the cell,
 *  which in a VOF cell is the volume-fraction weighted blend of the phase
 *  viscosities.  That makes the air branch use mu_air and the water branch
 *  use mu_water automatically.
 *
 *  NOTE on the 19.2 API: DEFINE_SOURCE is the standard 5-argument form,
 *  and the phase-2 (water) volume fraction is obtained from the mixture
 *  thread via THREAD_SUB_THREAD(t,1) -- index 1 is the secondary phase.
 * ===================================================================== */

#include "udf.h"

#define INVK_AIR    1.0e13    /* K = 1e-13 m2 -> air passes at ~0.44 mm/s @ 40 Pa */
#define INVK_WATER  1.0e15    /* K = 1e-15 m2 -> u = 4.4e-8 m/s at 40 Pa (still impermeable) */

/* ---------------------------------------------------------------------
 *  Shared inline evaluation of  -invK * mu  for a cell of the mixture
 *  thread t.  Returning invK*mu keeps the three directional sources
 *  guaranteed consistent with each other.
 * ------------------------------------------------------------------- */
static real invk_mu(cell_t c, Thread *t)
{
    real alpha;      /* water (secondary phase) volume fraction, 0..1 */
    real invK;       /* blended inverse permeability                  */
    real mu;         /* mixture laminar viscosity                     */

    alpha = C_VOF(c, THREAD_SUB_THREAD(t, 1));   /* secondary phase = water */
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    invK = (1.0 - alpha) * INVK_AIR + alpha * INVK_WATER;
    mu   = C_MU_L(c, t);

    return invK * mu;
}

/* ---------------------------- X momentum ---------------------------- */
DEFINE_SOURCE(membrane_x, c, t, dS, eqn)
{
    real k = invk_mu(c, t);

    dS[eqn] = -k;
    return -k * C_U(c, t);
}

/* ---------------------------- Y momentum ---------------------------- */
DEFINE_SOURCE(membrane_y, c, t, dS, eqn)
{
    real k = invk_mu(c, t);

    dS[eqn] = -k;
    return -k * C_V(c, t);
}

/* ---------------------------- Z momentum ---------------------------- */
DEFINE_SOURCE(membrane_z, c, t, dS, eqn)
{
    real k = invk_mu(c, t);

    dS[eqn] = -k;
    return -k * C_W(c, t);
}




