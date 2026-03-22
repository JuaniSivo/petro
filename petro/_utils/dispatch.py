"""
petro/_utils/dispatch.py
========================
Core dispatch utilities — the foundation of petro's type system.

Every petro function that accepts UnitFloat | ProbUnitFloat inputs calls
these helpers to:
  1. detect whether any input is probabilistic           (_is_prob)
  2. find the reference type that governs output         (_most_uncertain)
  3. unit-convert any quantia type                       (_coerce)
  4. find the common sample count                        (_n_samples)
  5. extract a flat list of floats in a target unit      (_as_samples)
  6. package a result back into the right quantia type   (_wrap)

Usage pattern inside a correlation
-----------------------------------
def rs_standing(P_sp, T_sp, SG_g, SG_o):
    ref = _most_uncertain(P_sp, T_sp, SG_g, SG_o)
    n   = _n_samples(P_sp, T_sp, SG_g, SG_o)

    P_vals  = _as_samples(P_sp,  "psi", n)
    T_vals  = _as_samples(T_sp,  "°F",  n)
    SG_g_v  = _as_samples(SG_g,  "1",   n)
    SG_o_v  = _as_samples(SG_o,  "1",   n)

    results = [_formula(P, T, g, o)
               for P, T, g, o in zip(P_vals, T_vals, SG_g_v, SG_o_v)]

    return _wrap(results, "scf_res/STB", ref)

This produces a UnitFloat when all inputs are UnitFloat, and a
ProbUnitFloat when any input is ProbUnitFloat — with zero branching in
the calling code.
"""
from __future__ import annotations
import array as _array
from quantia import UnitFloat, ProbUnitFloat, UnitArray, ProbUnitArray


def _is_prob(*args) -> bool:
    """Return True if any argument is ProbUnitFloat or ProbUnitArray."""
    return any(isinstance(a, (ProbUnitFloat, ProbUnitArray)) for a in args)


def _most_uncertain(*args) -> UnitFloat | ProbUnitFloat:
    """
    Return the most uncertain argument.

    Priority: ProbUnitFloat | ProbUnitArray > UnitFloat | UnitArray.

    Used to decide the output type: the result matches the most uncertain
    input so that deterministic inputs don't accidentally lose uncertainty.

    Raises
    ------
    TypeError  if no quantia type is found among args.
    """
    for a in args:
        if isinstance(a, (ProbUnitFloat, ProbUnitArray)):
            return a
    for a in args:
        if isinstance(a, (UnitFloat, UnitArray)):
            return a
    raise TypeError(
        f"_most_uncertain: no quantia type found among "
        f"{[type(a).__name__ for a in args]}"
    )


def _coerce(
    x: UnitFloat | ProbUnitFloat,
    target_unit: str,
) -> UnitFloat | ProbUnitFloat:
    """
    Convert x to target_unit, preserving type (UnitFloat → UnitFloat,
    ProbUnitFloat → ProbUnitFloat).

    A thin wrapper over quantia's .to() that makes the intent explicit
    in correlation code.
    """
    return x.to(target_unit)


def _n_samples(*args) -> int:
    """
    Return the common sample count from probabilistic arguments.

    Returns 1 if all arguments are deterministic (convention: one
    "sample" = the exact value).

    Raises
    ------
    ValueError  if multiple ProbUnitFloat arguments have different n.
    """
    ns = {a._n for a in args if isinstance(a, (ProbUnitFloat, ProbUnitArray))}
    if not ns:
        return 1
    if len(ns) > 1:
        raise ValueError(
            f"Mismatched sample counts among probabilistic inputs: {ns}. "
            "All ProbUnitFloat inputs in the same calculation must have "
            "the same number of samples. Use qu.config(n_samples=...) to "
            "ensure consistency."
        )
    return ns.pop()


def _as_samples(
    x: UnitFloat | ProbUnitFloat,
    target_unit: str,
    n: int,
) -> list[float]:
    """
    Return a list of n float values from x, converted to target_unit.

    UnitFloat  → [x.to(target_unit).value] * n   (broadcast)
    ProbUnitFloat → list(x.to(target_unit)._samples)

    This is the workhorse for building the input arrays that feed
    into correlation formulas.
    """
    x_conv = x.to(target_unit)
    if isinstance(x_conv, ProbUnitFloat):
        return list(x_conv._samples)
    return [x_conv.value] * n


def _wrap(
    value: float | list[float] | _array.array,
    unit: str,
    reference: UnitFloat | ProbUnitFloat,
) -> UnitFloat | ProbUnitFloat:
    """
    Package a computed result into the correct quantia output type.

    Parameters
    ----------
    value     : a single float (deterministic result) or a list / array
                of floats, one per MC sample.
    unit      : the unit of the result.
    reference : from _most_uncertain().  Determines the output type.

    Rules
    -----
    reference is ProbUnitFloat | ProbUnitArray
        → output is ProbUnitFloat
        → if value is a scalar, it is broadcast to reference._n samples

    reference is UnitFloat | UnitArray
        → output is UnitFloat
        → value must be a scalar float
    """
    if isinstance(reference, (ProbUnitFloat, ProbUnitArray)):
        if isinstance(value, (int, float)):
            arr = _array.array('d', [float(value)] * reference._n)
        else:
            arr = _array.array('d', value)
        return ProbUnitFloat._from_raw(arr, unit)
    return UnitFloat(float(value), unit)


__all__ = [
    "_is_prob",
    "_most_uncertain",
    "_coerce",
    "_n_samples",
    "_as_samples",
    "_wrap",
]