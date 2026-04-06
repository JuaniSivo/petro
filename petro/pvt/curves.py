# """
# petro/pvt/curves.py
# ===================
# Curve functions — thin wrappers over scalar correlations.

# Each function accepts a pressure array as the first argument and returns
# a matching array.  All other arguments are scalars (UnitFloat or
# ProbUnitFloat).

# Output type rule
# ----------------
# P_array is UnitArray of exact pressures.  Output type is determined
# by the most uncertain scalar argument:
#   · all UnitFloat  → UnitArray
#   · any ProbUnitFloat → ProbUnitArray

# The sample alignment guarantee holds: if SG_o is ProbUnitFloat(n=3000),
# sample k of SG_o is used at *every* pressure point.  The resulting
# ProbUnitArray captures how the full curve shifts across realizations.
# """
# from __future__ import annotations
# from quantia import UnitArray, ProbUnitArray, UnitFloat, ProbUnitFloat
# from petro._utils.dispatch import _is_prob
# from petro.pvt.correlations import (
#     rs_standing, bo_standing, bg, z_hall_yarborough, z_papay,
#     pseudo_critical_standing,
#     oil_viscosity_dead, oil_viscosity_saturated, gas_viscosity,
# )


# def _collect(elements: list) -> UnitArray | ProbUnitArray:
#     """Pack a list of UnitFloat or ProbUnitFloat into the right container."""
#     if isinstance(elements[0], ProbUnitFloat):
#         return ProbUnitArray(elements)
#     return UnitArray([e.value for e in elements], elements[0].unit)


# def rs_standing_curve(
#     P:    UnitArray | ProbUnitArray,
#     T:    UnitFloat | ProbUnitFloat,
#     SG_g: UnitFloat | ProbUnitFloat,
#     SG_o: UnitFloat | ProbUnitFloat,
# ) -> UnitArray | ProbUnitArray:
#     """Rs at each pressure in P [scf/STB]."""
#     return rs_standing(P, Rsb, T, SG_g, SG_o, Pb)


# def bo_standing_curve(
#     Rs:       UnitArray | ProbUnitArray,
#     SG_g:     UnitFloat | ProbUnitFloat,
#     SG_o:     UnitFloat | ProbUnitFloat,
#     T:        UnitFloat | ProbUnitFloat,
# ) -> UnitArray | ProbUnitArray:
#     """
#     Bo at each Rs in Rs_array [RB/STB].

#     Typically called with the output of rs_standing_curve so that the
#     pressure dependence is carried through correctly.
#     """
#     return bo_standing(Rs, SG_g, SG_o, T)


# def bg_curve(
#     P:    UnitArray | ProbUnitArray,
#     T:    UnitFloat | ProbUnitFloat,
#     SG_g: UnitFloat | ProbUnitFloat,
#     z_correlation: str = "hall_yarborough",
# ) -> UnitArray | ProbUnitArray:
#     """
#     Bg at each pressure in P [RB/scf].

#     Z-factor is computed internally from pseudo-critical properties.
#     """
#     Tpc, Ppc = pseudo_critical_standing(SG_g)
    
#     Ppr = P.to_si() / Ppc
#     Tpr = T.to_si() / Tpc
    
#     if z_correlation == "hall_yarborough":
#         z = z_hall_yarborough(Ppr, Tpr)
#     else:
#         z = z_papay(Ppr, Tpr)
    
#     return bg(P, T, z)


# def oil_viscosity_saturated_curve(
#     Rs:       UnitArray | ProbUnitArray,
#     api:      UnitFloat | ProbUnitFloat,
#     T:        UnitFloat | ProbUnitFloat,
# ) -> UnitArray | ProbUnitArray:
#     """Saturated-oil viscosity [cP] at each Rs in Rs_array."""
#     mu_od = oil_viscosity_dead(api, T)
#     return oil_viscosity_saturated(mu_od, Rs)


# def gas_viscosity_curve(
#     P:    UnitArray,
#     T:    UnitFloat | ProbUnitFloat,
#     SG_g: UnitFloat | ProbUnitFloat,
#     z_correlation: str = "hall_yarborough",
# ) -> UnitArray | ProbUnitArray:
#     """Gas viscosity [cP] at each pressure in P."""
#     Tpc, Ppc = pseudo_critical_standing(SG_g)
    
#     Ppr = P.to_si() / Ppc
#     Tpr = T.to_si() / Tpc
    
#     if z_correlation == "hall_yarborough":
#         z = z_hall_yarborough(Ppr, Tpr)
#     else:
#         z = z_papay(Ppr, Tpr)

#     return gas_viscosity(SG_g, T, P, z)