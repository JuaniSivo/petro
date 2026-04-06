"""
petro/pvt/models.py
===================
PVTFluid dataclass — the single object that carries all fluid description
through a project.

Design
------
PVTFluid stores fluid parameters (exact or probabilistic) and exposes
Rs(P), Bo(P), Bg(P), and curve variants.  When a field such as `api` is
ProbUnitFloat, every derived quantity is also ProbUnitFloat.

Serialization uses quantia's to_dict() / from_dict() for each field.
Strings and defaults round-trip through plain JSON.

Calibration factors (Rs_scale, Rs_offset, Bo_scale, Bo_offset) are
applied as:
    adjusted = scale × raw + offset
where scale is dimensionless "1" and offset is in the same unit as raw.
At default (scale=1, offset=0) the calibration is a no-op.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from quantia import UnitFloat, ProbUnitFloat, UnitArray, ProbUnitArray
from quantia._io import from_dict as _q_from_dict
from petro.pvt.correlations import (
    api_to_sg, pseudo_critical_standing, pbubble_point, bo_standing_below_pb,
    rs_bubble_point_from_gor, rs_standing, bo_standing, bg,
    z_hall_yarborough, z_papay, oil_compressibility,
    oil_viscosity_dead, oil_viscosity_saturated, gas_viscosity,
)
# from petro.pvt.curves import (
#     rs_standing_curve, bo_standing_curve, bg_curve,
#     oil_viscosity_saturated_curve, gas_viscosity_curve,
# )


@dataclass
class PVTFluid:
    """
    Complete fluid description.

    Parameters
    ----------
    fluid_type   : "black_oil" | "gas_condensate" | "gas"
    gor_sep_i    : initial separator GOR [m3_g_sc/m3_o_sc | scf/STB]
    api          : API gravity  [°API]
    SG_g         : gas specific gravity  ["1"]
    T_res        : reservoir temperature  [°C | °F | K | °R]
    Pb           : bubble-point pressure  [psia | bara | kg/cm2 | Pa]
    Rsb          : bubble-point Rs  [m3_g_sc/m3_o_sc | scf/STB]
    Bob          : bubble-point Bo  [m3_o_res/m3_o_sc | RB/STB]
    rho_o_b      : bubble-point oil density at reservoir conditions [kg/m3 | lb/ft3]
    c_o_b        : bubble-point oil compresibility [1/Pa, 1/bara, 1/psia]
    correlation  : "standing"
    z_correlation: "hall_yarborough" | "papay"
    Pb_scale     : 
    Pb_offset    : 
    Rs_scale     : multiplicative calibration factor  ["1"]
    Rs_offset    : additive calibration offset  ["scf/STB"]
    Bo_scale     : multiplicative calibration factor  ["1"]
    Bo_offset    : additive calibration offset  ["RB/STB"]
    """

    fluid_type:   str
    gor_sep_i:    UnitFloat | ProbUnitFloat
    api:          UnitFloat | ProbUnitFloat
    SG_g:         UnitFloat | ProbUnitFloat
    T_res:        UnitFloat | ProbUnitFloat
    P_sep:        UnitFloat | ProbUnitFloat
    T_sep:        UnitFloat | ProbUnitFloat

    correlation:   str = "standing"
    z_correlation: str = "hall_yarborough"         # "papay" | "hall_yarborough"

    Pb_scale:    UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    Pb_offset:   UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "Pa"))
    Rs_scale:    UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    Rs_offset:   UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "m3_g_sc/m3_o_sc"))
    Bo_scale:    UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    Bo_offset:   UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "m3_o_res/m3_o_sc"))
    Bg_scale:    UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    Bg_offset:   UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "m3_g_res/m3_g_sc"))
    mu_o_scale:  UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    mu_o_offset: UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "cP"))
    mu_g_scale:  UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    mu_g_offset: UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "cP"))
    c_o_scale:   UnitFloat = field(default_factory=lambda: UnitFloat(1.0, "1"))
    c_o_offset:  UnitFloat = field(default_factory=lambda: UnitFloat(0.0, "1/Pa"))

    # ── Derived helpers ───────────────────────────────────────────────────────

    @property
    def _SG_o(self) -> UnitFloat | ProbUnitFloat:
        return api_to_sg(self.api)

    def _z_factor(
        self, P: UnitFloat | ProbUnitFloat
    ) -> UnitFloat | ProbUnitFloat:
        Tpc, Ppc = pseudo_critical_standing(self.SG_g)
        P_pr = P.to("Pa") / Ppc
        T_pr = self.T_res.to("K") / Tpc
        if self.z_correlation == "papay":
            return z_papay(P_pr, T_pr)
        return z_hall_yarborough(P_pr, T_pr)

    @property
    def _Pb(self) -> UnitFloat | ProbUnitFloat:
        """Pressure at bubble point, calibration applied."""
        if self.correlation == "standing":
            raw = pbubble_point(self._Rsb, self.T_res, self.SG_g, self._SG_o)
        else:
            raise ValueError(f"Unknown correlation: '{self.correlation}'")
        return raw.to_si() * self.Pb_scale + self.Pb_offset

    @property
    def _Rsb(self) -> UnitFloat | ProbUnitFloat:
        """Pressure at bubble point, calibration applied."""
        if self.correlation == "standing":
            raw = rs_bubble_point_from_gor(self.gor_sep_i, self._SG_o, self. SG_g, self.P_sep, self.T_sep, include_stock_tank_gor=True)
        else:
            raise ValueError(f"Unknown correlation: '{self.correlation}'")
        return raw.to("m3_g_sc/m3_o_sc")

    @property
    def _Bob(self) -> UnitFloat | ProbUnitFloat:
        """Pressure at bubble point, calibration applied."""
        if self.correlation == "standing":
            raw = bo_standing_below_pb(self._Rsb, self.SG_g, self._SG_o, self.T_res)
        else:
            raise ValueError(f"Unknown correlation: '{self.correlation}'")
        return raw.to("m3_o_res/m3_o_sc")

    # ── Scalar methods ───────────────────────────────────────────────────────

    def Rs(self, P: UnitFloat | ProbUnitFloat) -> UnitFloat | ProbUnitFloat:
        """Solution GOR at pressure P [scf/STB], calibration applied."""
        if self.correlation == "standing":
            raw = rs_standing(P, self._Rsb, self.T_res, self.SG_g, self._SG_o, self._Pb)
        else:
            raise ValueError(f"Unknown correlation: '{self.correlation}'")
        return raw * self.Rs_scale.value + self.Rs_offset

    def Bo(self, P: UnitFloat | ProbUnitFloat) -> UnitFloat | ProbUnitFloat:
        """Oil FVF at pressure P [RB/STB], calibration applied."""
        if self.correlation == "standing":
            raw = bo_standing(P, self._Pb, self.Rs(P), self._Bob, self.SG_g, self._SG_o, self.T_res, self.c_o(P))
        else:
            raise ValueError(f"Unknown correlation: '{self.correlation}'")
        return raw * self.Bo_scale.value + self.Bo_offset

    def Bg(self, P: UnitFloat | ProbUnitFloat) -> UnitFloat | ProbUnitFloat:
        """Gas FVF at pressure P [RB/scf]."""
        z = self._z_factor(P)
        return bg(P, self.T_res, z)

    def mu_o(self, P: UnitFloat | ProbUnitFloat) -> UnitFloat | ProbUnitFloat:
        """Saturated oil viscosity at pressure P [cP]."""
        mu_dead = oil_viscosity_dead(self._SG_o, self.T_res)
        return oil_viscosity_saturated(mu_dead, self.Rs(P))

    def mu_g(self, P: UnitFloat | ProbUnitFloat) -> UnitFloat | ProbUnitFloat:
        """Gas viscosity at pressure P [cP]."""
        z = self._z_factor(P)
        return gas_viscosity(self.SG_g, self.T_res, P, z)

    def c_o(self, P: UnitFloat | ProbUnitFloat) -> UnitFloat | ProbUnitFloat:
        """Solution GOR at pressure P [scf/STB], calibration applied."""
        if self.correlation == "standing":
            raw = oil_compressibility(self.Rs(P), self.T_res, self.SG_g, self._SG_o, P, self._Pb)
        else:
            raise ValueError(f"Unknown correlation: '{self.correlation}'")
        return raw * self.c_o_scale.value + self.c_o_offset
    

    # ── Curve methods ─────────────────────────────────────────────────────────

    def Rs_curve(self, P_array: UnitArray) -> UnitArray | ProbUnitArray:
        """Rs at each pressure in P_array [m3/m3]."""
        return self.Rs(P_array)

    def Bo_curve(self, P_array: UnitArray) -> UnitArray | ProbUnitArray:
        """Bo at each pressure in P_array [m3/m3]."""
        return self.Bo(P_array)

    def Bg_curve(self, P_array: UnitArray) -> UnitArray | ProbUnitArray:
        """Bg at each pressure in P_array [m3/m3]."""
        return self.Bg(P_array)

    def mu_o_curve(self, P_array: UnitArray) -> UnitArray | ProbUnitArray:
        """Saturated oil viscosity at each pressure [cP]."""
        return self.mu_o(P_array)

    def mu_g_curve(self, P_array: UnitArray) -> UnitArray | ProbUnitArray:
        """Gas viscosity at each pressure [cP]."""
        return self.mu_g(P_array)

    def c_o_curve(self, P_array: UnitArray) -> UnitArray | ProbUnitArray:
        """Oil compressibility at each pressure [1/Pa]."""
        return self.c_o(P_array)

    # ── Serialization ─────────────────────────────────────────────────────────

    def to_dict(self) -> dict:
        def _s(x):
            return x.to_dict() if hasattr(x, "to_dict") else x
        return {
            "type":          "PVTFluid",
            "fluid_type":    self.fluid_type,
            "gor_sep_i":     _s(self.gor_sep_i),
            "api":           _s(self.api),
            "SG_g":          _s(self.SG_g),
            "T_res":         _s(self.T_res),
            "P_sep":         _s(self.P_sep),
            "T_sep":         _s(self.T_sep),
            "correlation":   self.correlation,
            "z_correlation": self.z_correlation,
            "Pb_scale":      _s(self.Pb_scale),
            "Pb_offset":     _s(self.Pb_offset),
            "Rs_scale":      _s(self.Rs_scale),
            "Rs_offset":     _s(self.Rs_offset),
            "Bo_scale":      _s(self.Bo_scale),
            "Bo_offset":     _s(self.Bo_offset),
            "Bg_scale":      _s(self.Bg_scale),
            "Bg_offset":     _s(self.Bg_offset),
            "mu_o_scale":    _s(self.mu_o_scale),
            "mu_o_offset":   _s(self.mu_o_offset),
            "mu_g_scale":    _s(self.mu_g_scale),
            "mu_g_offset":   _s(self.mu_g_offset),
            "c_o_scale":     _s(self.c_o_scale),
            "c_o_offset":    _s(self.c_o_offset),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "PVTFluid":
        if d.get("type") != "PVTFluid":
            raise ValueError(f"Expected type 'PVTFluid', got {d.get('type')!r}")
        return cls(
            fluid_type    = d["fluid_type"],
            gor_sep_i     = _q_from_dict(d["gor_sep_i"]),
            api           = _q_from_dict(d["api"]),
            SG_g          = _q_from_dict(d["SG_g"]),
            T_res         = _q_from_dict(d["T_res"]),
            P_sep         = _q_from_dict(d["P_sep"]),
            T_sep         = _q_from_dict(d["T_sep"]),
            correlation   = d.get("correlation",   "standing"),
            z_correlation = d.get("z_correlation", "hall_yarborough"),
            Pb_scale      = _q_from_dict(d.get("Pb_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            Pb_offset     = _q_from_dict(d.get("Pb_offset", {"type": "UnitFloat", "value": 1.0, "unit": "Pa"})),
            Rs_scale      = _q_from_dict(d.get("Rs_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            Rs_offset     = _q_from_dict(d.get("Rs_offset", {"type": "UnitFloat", "value": 1.0, "unit": "m3_g_sc/m3_o_sc"})),
            Bo_scale      = _q_from_dict(d.get("Bo_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            Bo_offset     = _q_from_dict(d.get("Bo_offset", {"type": "UnitFloat", "value": 1.0, "unit": "m3_o_res/m3_o_sc"})),
            Bg_scale      = _q_from_dict(d.get("Bg_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            Bg_offset     = _q_from_dict(d.get("Bg_offset", {"type": "UnitFloat", "value": 1.0, "unit": "m3_g_res/m3_g_sc"})),
            mu_o_scale    = _q_from_dict(d.get("mu_o_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            mu_o_offset   = _q_from_dict(d.get("mu_o_offset", {"type": "UnitFloat", "value": 1.0, "unit": "cP"})),
            mu_g_scale    = _q_from_dict(d.get("mu_g_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            mu_g_offset   = _q_from_dict(d.get("mu_g_offset", {"type": "UnitFloat", "value": 1.0, "unit": "cP"})),
            c_o_scale     = _q_from_dict(d.get("c_o_scale", {"type": "UnitFloat", "value": 1.0, "unit": "1"})),
            c_o_offset    = _q_from_dict(d.get("c_o_offset", {"type": "UnitFloat", "value": 1.0, "unit": "1/Pa"})),
        )