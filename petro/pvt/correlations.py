"""
petro/pvt/correlations.py
=========================
Standing (1947) PVT correlation family — scalar functions.

All functions:
  · Accept UnitFloat or ProbUnitFloat for any argument.
  · Return UnitFloat when all inputs are exact; ProbUnitFloat otherwise.
  · Convert units internally — callers pass any compatible unit.
  · Never do bare float arithmetic on UnitFloat objects.

Standing (1947) applicability ranges
-------------------------------------
Rs  :   100.0    -   2900.0   scf/STB
Bo  :     1.024  -      2.15  RB/STB
T   :   100.0    -    258.0   °F
API :    16.5    -     63.8   °API
SG_g:     0.59   -      0.95

Z-factor correlations
---------------------
Papay (1985)           — simple, explicit; good for low P_pr (< 7)
Hall-Yarborough (1974) — Newton-Raphson; industry standard for P_pr < 15

Gas viscosity
-------------
Lee-Gonzalez-Eakin (1966) — valid for 100 - 340°F, 200 - 8000 psia.

Temperature convention
-----------------------
All functions that express °F or °R in the published formula convert
internally via quantia's AffineUnit (.to("K") → multiply by 1.8 for °R,
or .to("°F") for explicit °F use). Callers may pass any temperature unit.
"""

from __future__ import annotations
from quantia import UnitFloat, UnitArray, ProbUnitFloat, ProbUnitArray, sg_to_api, api_to_sg
from quantia import math as qmath
from petro._utils.dispatch import (
    _most_uncertain, _n_samples, _as_samples, _wrap,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def pseudo_critical_standing(
    SG_g: UnitFloat | ProbUnitFloat
) -> tuple[UnitFloat | ProbUnitFloat, UnitFloat | ProbUnitFloat]:
    """
    Standing (1977) pseudo-critical properties for natural gas.

        Tpc [°R]   = 187 + 330  * SG_g - 71.5 * SG_g ** 2
        Ppc [psia] = 706 - 51.7 * SG_g - 11.1 * SG_g ** 2

    Returns
    -------
    (Tpc, Ppc) : both UnitFloat | ProbUnitFloat
    Tpc in "K", Ppc in "Pa".
    """

    Tpc = (187 + 330  * SG_g - 71.5 * SG_g * SG_g) * UnitFloat(1, "°R")
    Ppc = (706 - 51.7 * SG_g - 11.1 * SG_g * SG_g) * UnitFloat(1, "psia")

    return Tpc.to_si(), Ppc.to_si()


def stock_tank_gor(
    SGo: UnitFloat | ProbUnitFloat,
    SGg_sep: UnitFloat | ProbUnitFloat,
    p_sp: UnitFloat | ProbUnitFloat,
    T_sp: UnitFloat | ProbUnitFloat
) -> UnitFloat | ProbUnitFloat:
    # Stock-tank GOR
    a1 =  0.3818
    a2 = -5.506
    a3 =  2.902
    a4 =  1.327
    a5 = -0.7355

    log_gor_st = a1 \
               + a2 * qmath.log(SGo) \
               + a3 * qmath.log(SGg_sep) \
               + a4 * qmath.log(p_sp.to("psia")) \
               + a5 * qmath.log(T_sp.to("°F"))
    
    gor_st = (10 ** log_gor_st) * UnitFloat(1, "scf/STB")

    return gor_st.to_si() * UnitFloat(1, "m3_g_sc/m3_o_sc")


def rs_bubble_point_from_gor_sep(gor_sep):
    # separator GOR
    return gor_sep.to_si() * UnitFloat(1, "m3_g_sc/m3_o_sc")


def rs_bubble_point_from_gor(gor_sep, SGo=None, SGg_sep=None, p_sp=None, T_sp=None, include_stock_tank_gor=False):
    rs = rs_bubble_point_from_gor_sep(gor_sep)
    if include_stock_tank_gor:
        rs  = rs + stock_tank_gor(SGo, SGg_sep, p_sp, T_sp)
    
    return rs.to("m3_g_sc/m3_o_sc")


def rs_bubble_point_from_pb(Pb, T, SG_g, SG_o):
    T_adim = T.to("°F") * UnitFloat(1, "1/°F")
    api_adim = sg_to_api(SG_o) * UnitFloat(1, "1/°API")
    P_adim = Pb.to("psia") * UnitFloat(1, "1/psia")
    
    exponent = 0.00091 * T_adim - 0.0125 * api_adim
    Cpb  = P_adim / 18.2 + 1.4
    Rsb =  SG_g * (Cpb / (10 ** exponent)) ** (1 / 0.83) * UnitFloat(1, "scf/STB")

    return Rsb.to("m3_g_sc/m3_o_sc")


def pbubble_point(Rsb, T, SG_g, SG_o):
    T_adim = T.to("°F") * UnitFloat(1, "1/°F")
    api_adim = sg_to_api(SG_o) * UnitFloat(1, "1/°API")
    Rsb_adim = Rsb.to("scf/STB") * UnitFloat(1, "STB/scf")

    exponent = 0.00091 * T_adim - 0.0125 * api_adim
    Cpb = (Rsb_adim / SG_g) ** 0.83 * 10 ** exponent
    Pb = 18.2 * (Cpb - 1.4) * UnitFloat(1, "psia")

    return Pb.to_si()


# ─────────────────────────────────────────────────────────────────────────────
# Solution GOR — rs_standing
# ─────────────────────────────────────────────────────────────────────────────

def rs_standing(P, Rsb, T, SG_g, SG_o, Pb):
    """
    Standing (1947) solution GOR correlation.

        Rs = SG_g x [(P/18.2 + 1.4) x 10^(0.0125·API - 0.00091·T_F)]^1.205

    Parameters
    ----------
    P    : pressure  [any pressure unit → converted to psi]
    T    : temperature  [any temperature unit → converted to °F]
    SG_g : gas specific gravity  ["1"]
    SG_o : oil specific gravity  ["1"]  — API derived internally

    Returns
    -------
    Rs : UnitFloat | ProbUnitFloat  [m3/m3 = adimensional]
    """
    T_adim = T.to("°F") * UnitFloat(1, "1/°F")
    api_adim = sg_to_api(SG_o) * UnitFloat(1, "1/°API")
    P_adim = P.to("psia") * UnitFloat(1, "1/psia")
    Pb_adim = Pb.to("psia") * UnitFloat(1, "1/psia")
    
    exponent = 0.00091 * T_adim - 0.0125 * api_adim
    Cp  = P_adim / 18.2 + 1.4
    Rs_below_pb =  SG_g * (Cp / (10 ** exponent)) ** (1 / 0.83) * UnitFloat(1, "scf/STB")
    


    mask = 0.5 * ((P-Pb)/abs(P-Pb) + 1)
    Rs = (1-mask) * Rs_below_pb.to("m3_g_sc/m3_o_sc") \
         + mask   * Rsb

    return Rs


# ─────────────────────────────────────────────────────────────────────────────
# Oil FVF — bo_standing
# ─────────────────────────────────────────────────────────────────────────────

def bo_standing_below_pb(Rs, SG_g, SG_o, T):
    Rs_adim = Rs.to("scf/STB") * UnitFloat(1, "STB/scf")
    T_adim = T.to("°F") * UnitFloat(1, "1/°F")

    Cbo = Rs_adim * (SG_g / SG_o) ** 0.5 + 1.25 * T_adim
    Bo = (0.9759 + 12 * 10e-6 * (Cbo ** 1.2)) * UnitFloat(1, "RB/STB")
    
    return Bo


def bo_standing_above_pb(P, Pb, c_o, Bob):
    Bo = Bob * qmath.exp(c_o.to("1/psia") * (Pb.to("psia") - P.to("psia")))
    
    return Bo


def bo_standing(P, Pb, Rs, Bob, SG_g, SG_o, T, c_o):
    """
    Standing (1947) oil formation volume factor.

        F  = Rs x √(SG_g / SG_o) + 1.25 x T_F
        Bo = 0.972 + 1.47 x 10⁻⁴ x F^1.175

        Bob = 0.9759 + 12 * 10e-5 * CBob ** 1.2
        Cbob = Rs x (SGg/SGo) ** 0.5 + 1.25 x T

    Parameters
    ----------
    Rs   : solution GOR  [scf/STB]
    SG_g : gas specific gravity  ["1"]
    SG_o : oil specific gravity  ["1"]
    T    : temperature  [any temperature unit → converted to °F]

    Returns
    -------
    Bo : UnitFloat | ProbUnitFloat  [RB/STB]
    """
    mask = 0.5 * ((P-Pb)/abs(P-Pb) + 1)
    Bo = (1-mask) * bo_standing_below_pb(Rs, SG_g, SG_o, T) \
         + mask   * bo_standing_above_pb(P, Pb, c_o, Bob)
    
    return Bo


# ─────────────────────────────────────────────────────────────────────────────
# Gas Z-factor
# ─────────────────────────────────────────────────────────────────────────────

def z_papay(P_pr, T_pr):
    """Papay (1985) explicit z-factor. Valid for P_pr < ~7."""
    return (1.0
            - 3.52  * P_pr / (10.0 ** (0.9813 * T_pr))
            + 0.274 * P_pr * P_pr / (10.0 ** (0.8157 * T_pr)))


def z_hall_yarborough(P_pr, T_pr):
    """
    Hall-Yarborough (1974) z-factor via Newton-Raphson.

    Solves F(y) = 0 for reduced density y, then z = alpha·P_pr / y.

    Variables
    ---------
    t = 1 / T_pr
    alpha = 0.06125·t·exp(-1.2·(1-t)²)          [positive]
    b = 14.76t - 9.76t² + 4.58t³
    c = 90.7t - 242.2t² + 42.4t³
    d = 2.18 + 2.82t

    F(y) = -alpha·Ppr + (y+y²+y³-y⁴)/(1-y)³ - b·y² + c·y^d

    F'(y) = (1+4y+4y²-4y³+y⁴)/(1-y)⁴ - 2b·y + c·d·y^(d-1)

    Valid range: 1.2 ≤ T_pr, 0.2 ≤ P_pr ≤ 15.
    """

    t = 1.0 / T_pr
    alpha = 0.06125 * t    * qmath.exp(-1.2 * (1.0 - t) ** 2)
    b     = 14.76   * t    - 9.76  * t * t + 4.58  * t ** 3
    c     = 90.7    * t    - 242.2 * t * t + 42.4  * t ** 3
    d     = 2.18    + 2.82 * t

    # Initial guess: first-order linear approximation F(y) ≈ -α·Ppr + y = 0
    y = alpha * P_pr

    for _ in range(60):
        one_minus_y    = 1.0 - y
        omy3           = one_minus_y ** 3
        omy4           = one_minus_y * omy3

        F  = (- alpha * P_pr
              + (y + y*y + y**3 - y**4) / omy3
              - b * y * y
              + c * y ** d)

        dF = ((1.0 + 4*y + 4*y*y - 4*y**3 + y**4) / omy4
              - 2.0 * b * y
              + c * d * y ** (d - 1.0))

        step  = F / dF
        y_new = y - step

        if isinstance(y_new, UnitFloat):
            # guard against negative density
            if y_new <= 0.0: mask = UnitFloat(1, "1")
            if y_new >  0.0: mask = UnitFloat(0, "1")
            y = (1 - mask) * y_new + mask * y * 0.5
            if step/y < 1e-9:
                break
        if isinstance(y_new, UnitArray):
            mask = UnitArray([1 if i<=0 else 0 for i in y_new._data], "1")
            y = (1 - mask) * y_new + mask * y * 0.5
            if (step/y).max() < 1e-9:
                break
        if isinstance(y_new, ProbUnitFloat):
            mask = ProbUnitFloat._from_raw([1 if i<=0 else 0 for i in y_new._samples], "1")
            y = (1 - mask) * y_new + mask * y * 0.5
            if (step/y).max() < 1e-9:
                break
        if isinstance(y_new, ProbUnitArray):
            mask = ProbUnitArray._from_flat([1 if i<=0 else 0 for i in y_new._data], "1", y_new._len, y_new._n)
            y = (1 - mask) * y_new + mask * y * 0.5
            if max([puf.max() for puf in (step/y)]) < 1e-9:
                break

    return alpha * P_pr / y


# ─────────────────────────────────────────────────────────────────────────────
# Gas FVF — bg
# ─────────────────────────────────────────────────────────────────────────────

def bg(P, T, z):
    """
    Gas formation volume factor.

    Derived from real gas law at standard conditions
    (P_sc = 14.696 psia, T_sc = 60 °F = 519.67 °R):

        Bg [RB/scf] = (P_sc / (T_sc x 5.615)) x z x T_R / P
                    = 0.005035 x z x T_R / P

    where T_R = T [°R] = T [K] x 9/5.

    Parameters
    ----------
    P : pressure  [any pressure unit → converted to psi]
    T : temperature  [any temperature unit → converted to K → °R]
    z : z-factor  ["1"]

    Returns
    -------
    Bg : UnitFloat | ProbUnitFloat  [RB/scf]
    """

    T_sc = UnitFloat(60, "°F")
    p_sc = UnitFloat(1, "atm")

    Bg = z * p_sc.to_si() * T.to_si() / (T_sc.to_si() * P.to_si())

    return Bg.to("m3_g_res/m3_g_sc")


# ─────────────────────────────────────────────────────────────────────────────
# Oil Viscosity
# ─────────────────────────────────────────────────────────────────────────────

def oil_viscosity_dead(SG_o, T):
    """
    Beal (1946) dead-oil viscosity correlation.

        A    = 10^(0.43 + 8.33/API)
        μ_od = (0.32 + 1.8x10⁷ / API^4.53) x (360 / (T_F + 200))^A

    Parameters
    ----------
    api : API gravity  [°API]
    T   : temperature  [any temperature unit → converted to °F]

    Returns
    -------
    μ_od : UnitFloat | ProbUnitFloat  [cP]
    """

    api_adim = sg_to_api(SG_o) * UnitFloat(1, "1/°API")
    T_adim = T.to("°F") * UnitFloat(1, "1/°F")

    A     = 10.0 ** (0.43 + 8.33 / api_adim)
    mu_od = (0.32 + 1.8e7 / (api_adim ** 4.53)) * (360.0 / (T_adim + 200.0)) ** A
    mu_od = mu_od * UnitFloat(1, "cP")

    return mu_od


def oil_viscosity_saturated(mu_od, Rs):
    """
    Chew-Connally (1959) saturated-oil viscosity correlation.

        A   = 10^[Rs · (2.2x10⁻⁷·Rs - 7.4x10⁻⁴)]
        b   = 0.68/10^(8.62x10⁻⁵·Rs)
            + 0.25/10^(1.1x10⁻³·Rs)
            + 0.062/10^(1.082x10⁻²·Rs)
        μ_o = A · μ_od^b

    Parameters
    ----------
    mu_od : dead-oil viscosity  [cP]
    Rs    : solution GOR  [scf/STB]

    Returns
    -------
    μ_o : UnitFloat | ProbUnitFloat  [cP]
    """

    Rs_adim = Rs.to("scf/STB") * UnitFloat(1, "STB/scf")
    mu_od_adim = mu_od.to("cP") * UnitFloat(1, "1/cP")
    
    log_A = Rs_adim * (2.2e-7 * Rs_adim - 7.4e-4)
    A     = 10.0 ** log_A
    b     = (0.68  / 10.0 ** (8.62e-5   * Rs_adim)
           + 0.25  / 10.0 ** (1.1e-3    * Rs_adim)
           + 0.062 / 10.0 ** (1.082e-2  * Rs_adim))
    mu_o  = A * mu_od_adim ** b
    mu_o  = mu_o * UnitFloat(1, "cP")
    
    return mu_o


# ─────────────────────────────────────────────────────────────────────────────
# Gas Viscosity
# ─────────────────────────────────────────────────────────────────────────────

def gas_viscosity(SG_g, T, P, z):
    """
    Lee-Gonzalez-Eakin (1966) gas viscosity correlation.

        M_g  = 28.97 · SG_g
        ρ_g  = P · M_g / (z · 10.73 · T_R)   [lb/ft³]
               x 0.016018                      [g/cc]
        K    = (9.4 + 0.02·M_g) · T_R^1.5 / (209 + 19·M_g + T_R)
        X    = 3.5 + 986/T_R + 0.01·M_g
        Y    = 2.4 - 0.2·X
        μ_g  = K · 10⁻⁴ · exp(X · ρ_g^Y)   [cP]

    Valid range: 100-340 °F, 200-8 000 psia.

    Parameters
    ----------
    SG_g : gas specific gravity  ["1"]
    T    : temperature  [any temperature unit → converted to K → °R]
    P    : pressure  [any pressure unit → converted to psi]
    z    : z-factor  ["1"]

    Returns
    -------
    μ_g : UnitFloat | ProbUnitFloat  [cP]
    """
    
    P_adim = P.to("psia") * UnitFloat(1, "1/psia")
    T_adim = T.to("°R") * UnitFloat(1, "1/°R")
    
    m_g   = 28.97 * SG_g
    rho_g = P_adim * m_g / (z * 10.73 * T_adim) * 0.016018
    K     = (9.4 + 0.02 * m_g) * T_adim ** 1.5 / (209.0 + 19.0 * m_g + T_adim)
    X     = 3.5 + 986.0 / T_adim + 0.01 * m_g
    Y     = 2.4 - 0.2 * X
    mu_g  = K * 1e-4 * qmath.exp(X * rho_g ** Y) * UnitFloat(1, "cP")

    return mu_g


# ─────────────────────────────────────────────────────────────────────────────
# Oil Density at Reservoir Conditions
# ─────────────────────────────────────────────────────────────────────────────

def oil_density_below_pb(Rs, Bo, SG_g, SG_o):
    rho_o_sc = SG_o * UnitFloat(1000 , "kg/m^3")

    rho_o_sc_adim = rho_o_sc.to("lb/ft3") * UnitFloat(1, "ft3/lb")
    Rs_adim = Rs.to("scf/STB") * UnitFloat(1, "STB/scf")
    Bo_adim = Bo.to("RB/STB") * UnitFloat(1, "STB/RB")

    rho_o = UnitFloat(1, "lb/ft3") * (rho_o_sc_adim + 0.01357 * Rs_adim  * SG_g) / Bo_adim

    return rho_o.to("kg/m^3")


def oil_density_above_pb(rho_o_b, P, Pb, c_o):
    P_adim = P.to("psia") * UnitFloat(1, "1/psia")
    Pb_adim = Pb.to("psia") * UnitFloat(1, "1.psia")
    c_o_adim = c_o.to("1/psia") * UnitFloat(1, "psia")
    rho_o_b = rho_o_b.to("kg/m^3") * UnitFloat(1, "m^3/kg")

    rho_o = rho_o_b * qmath.exp(c_o_adim * (P_adim - Pb_adim))

    return rho_o


def oil_density(Rs, Bo, rho_o_b, SG_g, SG_o, P, Pb, c_o):
    """
    Below bubble point:
        rho_o = (rho_o_ST + 0.01357 * Rs  * SG_g) / Bo
    Above bubble point:
        rho_o = rho_o_b * exp(c_o * (P - Pb))

    Parameters
    ----------
    Rs   : 
    Rsb  : bubblepoint Rs
    Bo   : 
    Bob  : bubblepoint Bo
    P    : pressure [Any and convert to psia]
    Pb   : bubblepoint pressure [Any and convert to psia]
    SG_g : gas specific gravity  ["1"]
    SG_o : gas specific gravity  ["1" and convert to rho oil in lb/ft^3]
    c_o  : 

    Returns
    -------
    rho_o : UnitFloat | ProbUnitFloat  [kg/m^3]
    """
    mask = 0.5 * ((P-Pb)/abs(P-Pb) + 1)
    rho_o = (
            (1-mask) * oil_density_below_pb(Rs, Bo, SG_g, SG_o) \
            + mask   * oil_density_above_pb(rho_o_b, P, Pb, c_o)
            ) * UnitFloat(1, "lb/ft3")

    return rho_o.to("kg/m^3")


# ─────────────────────────────────────────────────────────────────────────────
# Coetficient or Isothermal Compressibility of Oil
# ─────────────────────────────────────────────────────────────────────────────

def oil_compressibility_below_pb(Rs, T, SG_o, P):
    Rs_adim  = Rs.to("scf/STB") * UnitFloat(1, "STB/scf")
    T_adim   = T.to("°F")       * UnitFloat(1, "1/°F")
    api_adim = sg_to_api(SG_o)  * UnitFloat(1, "1/°API")
    P_adim   = P.to("psia")     * UnitFloat(1, "1/psia")
    
    exponent = - 7.633 \
               - 1.497 * qmath.log(P_adim) \
               + 1.115 * qmath.log(T_adim) \
               + 0.533 * qmath.log(api_adim) \
               + 0.184 * qmath.log(Rs_adim)
    
    return qmath.exp(exponent) * UnitFloat(1, "1/psia")

def oil_compressibility_above_pb(Rs, T, SG_g, SG_o, P):
    
    Rs_adim  = Rs.to("scf/STB") * UnitFloat(1, "STB/scf")
    T_adim   = T.to("°F")       * UnitFloat(1, "1/°F")
    api_adim = sg_to_api(SG_o)  * UnitFloat(1, "1/°API")
    P_adim   = P.to("psia")     * UnitFloat(1, "1/psia")
    c_o_min  = UnitFloat(1e-7, "1/psia")

    num = - 1433 \
          + 5     * Rs_adim \
          + 17.2  * T_adim \
          + -1180 * SG_g \
          + 12.61 * api_adim
    
    c_o = num / (10e5 * P_adim) * UnitFloat(1, "1/psia")

    mask = 0.5 * ((c_o-c_o_min)/abs(c_o-c_o_min) + 1)
    c_o_adj = ((1-mask) * c_o_min + mask * c_o)

    return c_o_adj.to("1/Pa")

def oil_compressibility(Rs, T, SG_g, SG_o, P, Pb):
    """
    Lee-Gonzalez-Eakin (1966) gas viscosity correlation.

        M_g  = 28.97 · SG_g
        ρ_g  = P · M_g / (z · 10.73 · T_R)   [lb/ft³]
               x 0.016018                      [g/cc]
        K    = (9.4 + 0.02·M_g) · T_R^1.5 / (209 + 19·M_g + T_R)
        X    = 3.5 + 986/T_R + 0.01·M_g
        Y    = 2.4 - 0.2·X
        μ_g  = K · 10⁻⁴ · exp(X · ρ_g^Y)   [cP]

    Valid range: 100-340 °F, 200-8 000 psia.

    Parameters
    ----------
    SG_g : gas specific gravity  ["1"]
    T    : temperature  [any temperature unit → converted to K → °R]
    P    : pressure  [any pressure unit → converted to psi]
    z    : z-factor  ["1"]

    Returns
    -------
    μ_g : UnitFloat | ProbUnitFloat  [cP]
    """    
    mask = 0.5 * ((P-Pb)/abs(P-Pb) + 1)
    c_o = (
          (1-mask) * oil_compressibility_below_pb(Rs, T, SG_o, P) \
          + mask   * oil_compressibility_above_pb(Rs, T, SG_g, SG_o, P)
          )

    return c_o