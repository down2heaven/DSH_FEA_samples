/* =====================================================================
 *  membrane.c  -  Phase-selective Darcy momentum sink for the
 *                 breathable-waterproof membrane case (rev. 2026-11-23)
 *
 *  Target : ANSYS Fluent 19.2, 3ddp, VOF (2 phases: air = primary,
 *           water-liquid = secondary)
 *  Attach : cell zone "porous" -> mixture -> source terms,
 *           X/Y/Z momentum = membrane_x / membrane_y / membrane_z
 *
 *  Sink (Darcy, Forchheimer term neglected, pore Re << 1):
 *
 *      S_i = -R(alpha) * u_i            [N/m3]
 *      R   = (1-w) * MU_AIR/K_AIR + w * MU_WATER/K_WATER
 *      w   = 1 - exp(-alpha/ALPHA_C)
 *
 *  Why w and not alpha:
 *  The old version blended invK linearly in alpha and multiplied by the
 *  mixture viscosity.  A cell holding 1-10 % water then had roughly the
 *  air resistance, so the water front advanced into the membrane at an
 *  "air" rate and only stopped once the cell was full (observed in
 *  film13: front ~0.25 mm deep after 0.46 s, Darcy predicts ~0.05 um).
 *  With ALPHA_C = 5e-3 the resistance is already 86 % of the water value
 *  at alpha = 0.01 and 98 % at alpha = 0.02, so any cell that starts to
 *  wet is closed off immediately, while a dry cell keeps the air value.
 *
 *  Phase viscosities are taken as constants equal to the material
 *  properties in the case (air 1.789e-5, water-liquid 1.003e-3 Pa s);
 *  using the mixture viscosity would again mix the two branches.
 *
 *  Values (resistance R = mu/K):
 *    air   : K = 1e-13 m2 -> R = 1.79e8  Pa s/m2  (breathable)
 *    water : K = 1e-15 m2 -> R = 1.00e12 Pa s/m2  (u ~ 1e-7 m/s at 50 Pa
 *                                                  across 0.5 mm)
 *
 *  dS/du_i = -R is exact, so the sink is fully implicit in u_i.
 *  alpha is lagged (taken from the current iterate), which is the usual
 *  explicit treatment of property dependence.
 * ===================================================================== */

#include "udf.h"

#define K_AIR      1.0e-13    /* m2 */
#define K_WATER    1.0e-15    /* m2 */
#define MU_AIR     1.789e-5   /* Pa s, must match material "air"          */
#define MU_WATER   1.003e-3   /* Pa s, must match material "water-liquid" */
#define ALPHA_C    5.0e-3     /* wetting scale of the switch              */

static real resistance(cell_t c, Thread *t)
{
    real alpha = C_VOF(c, THREAD_SUB_THREAD(t, 1));   /* water fraction */
    real w;

    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    w = 1.0 - exp(-alpha / ALPHA_C);

    return (1.0 - w) * (MU_AIR / K_AIR) + w * (MU_WATER / K_WATER);
}

DEFINE_SOURCE(membrane_x, c, t, dS, eqn)
{
    real r = resistance(c, t);
    dS[eqn] = -r;
    return -r * C_U(c, t);
}

DEFINE_SOURCE(membrane_y, c, t, dS, eqn)
{
    real r = resistance(c, t);
    dS[eqn] = -r;
    return -r * C_V(c, t);
}

DEFINE_SOURCE(membrane_z, c, t, dS, eqn)
{
    real r = resistance(c, t);
    dS[eqn] = -r;
    return -r * C_W(c, t);
}
