"""Tests for petro._utils.dispatch — the type-dispatch foundation."""
import pytest
import quantia as qu
from petro._utils.dispatch import (
    _is_prob, _most_uncertain, _coerce,
    _n_samples, _as_samples, _wrap,
)


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def uf_m():
    return qu.Q(1000.0, "m")

@pytest.fixture
def uf_s():
    return qu.Q(10.0, "s")

@pytest.fixture
def puf_m():
    with qu.config(seed=0, n_samples=200):
        return qu.ProbUnitFloat.uniform(800.0, 1200.0, "m")

@pytest.fixture
def puf_s():
    with qu.config(seed=1, n_samples=200):
        return qu.ProbUnitFloat.uniform(8.0, 12.0, "s")


# ── _is_prob ──────────────────────────────────────────────────────────────────

def test_is_prob_all_exact(uf_m, uf_s):
    assert not _is_prob(uf_m, uf_s)

def test_is_prob_one_prob(uf_m, puf_s):
    assert _is_prob(uf_m, puf_s)

def test_is_prob_empty():
    assert not _is_prob()

def test_is_prob_non_quantia():
    assert not _is_prob(1.0, "hello")


# ── _most_uncertain ───────────────────────────────────────────────────────────

def test_most_uncertain_picks_prob(uf_m, puf_s):
    result = _most_uncertain(uf_m, puf_s)
    assert result is puf_s

def test_most_uncertain_all_exact(uf_m, uf_s):
    result = _most_uncertain(uf_m, uf_s)
    assert isinstance(result, qu.UnitFloat)

def test_most_uncertain_no_quantia():
    with pytest.raises(TypeError, match="no quantia type"):
        _most_uncertain(1.0, "hello", [1, 2, 3])


# ── _coerce ───────────────────────────────────────────────────────────────────

def test_coerce_unitfloat_converts(uf_m):
    result = _coerce(uf_m, "km")
    assert isinstance(result, qu.UnitFloat)
    assert result.value == pytest.approx(1.0)

def test_coerce_probunitfloat_converts(puf_m):
    result = _coerce(puf_m, "km")
    assert isinstance(result, qu.ProbUnitFloat)
    assert result.mean().value == pytest.approx(1.0, rel=0.05)
    assert result._n == puf_m._n

def test_coerce_incompatible_raises(uf_m):
    from quantia import IncompatibleUnitsError
    with pytest.raises(IncompatibleUnitsError):
        _coerce(uf_m, "s")


# ── _n_samples ────────────────────────────────────────────────────────────────

def test_n_samples_all_exact(uf_m, uf_s):
    assert _n_samples(uf_m, uf_s) == 1

def test_n_samples_with_prob(uf_m, puf_m):
    assert _n_samples(uf_m, puf_m) == 200

def test_n_samples_mismatch_raises():
    with qu.config(seed=0):
        p1 = qu.ProbUnitFloat.uniform(0, 1, "1", n=100)
        p2 = qu.ProbUnitFloat.uniform(0, 1, "1", n=300)
    with pytest.raises(ValueError, match="Mismatched"):
        _n_samples(p1, p2)


# ── _as_samples ───────────────────────────────────────────────────────────────

def test_as_samples_unitfloat_broadcasts(uf_m):
    samples = _as_samples(uf_m, "km", 5)
    assert samples == [1.0, 1.0, 1.0, 1.0, 1.0]

def test_as_samples_unitfloat_converts():
    x = qu.Q(1.0, "km")
    samples = _as_samples(x, "m", 3)
    assert samples == [1000.0, 1000.0, 1000.0]

def test_as_samples_probunitfloat_returns_samples(puf_m):
    samples = _as_samples(puf_m, "m", 200)
    assert len(samples) == 200
    assert all(800.0 <= s <= 1200.0 for s in samples)

def test_as_samples_probunitfloat_converts_units(puf_m):
    samples = _as_samples(puf_m, "km", 200)
    assert all(0.8 <= s <= 1.2 for s in samples)


# ── _wrap ─────────────────────────────────────────────────────────────────────

def test_wrap_exact_reference_scalar(uf_m):
    result = _wrap(9.81, "m/s^2", uf_m)
    assert isinstance(result, qu.UnitFloat)
    assert result.value == pytest.approx(9.81)
    assert str(result.unit) == "m/s^2"

def test_wrap_prob_reference_broadcasts_scalar(puf_m):
    result = _wrap(9.81, "m/s^2", puf_m)
    assert isinstance(result, qu.ProbUnitFloat)
    assert result._n == 200
    assert all(v == pytest.approx(9.81) for v in result._samples)

def test_wrap_prob_reference_accepts_list(puf_m):
    samples = list(range(200))
    result = _wrap(samples, "m", puf_m)
    assert isinstance(result, qu.ProbUnitFloat)
    assert list(result._samples) == pytest.approx(samples)

def test_wrap_roundtrip_pattern(uf_m, puf_s):
    """Simulate the dispatch pattern used in every correlation."""
    ref = _most_uncertain(uf_m, puf_s)
    n   = _n_samples(uf_m, puf_s)
    d   = _as_samples(uf_m, "m", n)
    t   = _as_samples(puf_s, "s", n)
    speeds = [di / ti for di, ti in zip(d, t)]
    result = _wrap(speeds, "m/s", ref)
    assert isinstance(result, qu.ProbUnitFloat)
    assert result._n == 200
    assert result.mean().value == pytest.approx(100.0, rel=0.1)


# ── Temperature (AffineUnit) coercion ─────────────────────────────────────────
# Standing uses °F internally. Reservoir temperature is often entered in °C.
# _coerce must handle the affine path via quantia's .to().

def test_coerce_celsius_to_fahrenheit():
    t_c = qu.Q(100.0, "°C")
    t_f = _coerce(t_c, "°F")
    assert isinstance(t_f, qu.UnitFloat)
    assert t_f.value == pytest.approx(212.0, rel=1e-5)

def test_coerce_fahrenheit_to_kelvin():
    t_f = qu.Q(32.0, "°F")
    t_k = _coerce(t_f, "K")
    assert t_k.value == pytest.approx(273.15, rel=1e-5)

def test_coerce_prob_celsius_to_fahrenheit():
    with qu.config(seed=0, n_samples=200):
        t_c = qu.ProbUnitFloat.uniform(80.0, 120.0, "°C")
    t_f = _coerce(t_c, "°F")
    assert isinstance(t_f, qu.ProbUnitFloat)
    assert t_f._n == 200
    # All samples should be in [176, 248] °F
    assert all(176.0 <= s <= 248.0 for s in t_f._samples)

def test_as_samples_celsius_to_fahrenheit():
    t_c = qu.Q(100.0, "°C")
    vals = _as_samples(t_c, "°F", 3)
    assert vals == pytest.approx([212.0, 212.0, 212.0], rel=1e-5)

def test_as_samples_prob_celsius_to_fahrenheit():
    with qu.config(seed=0, n_samples=100):
        t_c = qu.ProbUnitFloat.uniform(80.0, 120.0, "°C")
    vals = _as_samples(t_c, "°F", 100)
    assert len(vals) == 100
    assert all(176.0 <= v <= 248.0 for v in vals)