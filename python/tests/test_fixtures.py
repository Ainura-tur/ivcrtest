"""Checks on the reference results recorded from R.

These tests need only numpy. They confirm that each fixture is complete and
internally consistent, and that the random draws R recorded reproduce what R
computed from them. The replay tests of the port will rely on that.
"""

import math
from statistics import NormalDist

import numpy as np
import pytest

from conftest import arr
from ivcrtest import R_REFERENCE_VERSION

FIXTURES = {
    "algebra": "exact",
    "nearpd": "exact",
    "data": "exact",
    "bootstrap": "replay",
    "union_bounds": "exact",
    "cihybrid": "replay",
    "end_to_end": "tolerance",
}


def cr_delta(rho):
    """delta = (rho_xz, h), h = (rho_xy rho_xz - rho_yz) / sqrt(1 - rho_xy^2)."""
    rho_xy, rho_xz, rho_yz = rho
    return np.array([rho_xz, (rho_xy * rho_xz - rho_yz) / math.sqrt(1 - rho_xy**2)])


def corr3(x, y, z):
    """(rho_xy, rho_xz, rho_yz), in the order the R package uses."""
    c = np.corrcoef(np.vstack([x, y, z]))
    return np.array([c[0, 1], c[0, 2], c[1, 2]])


@pytest.mark.parametrize("name, kind", FIXTURES.items())
def test_fixture_metadata(fixture, name, kind):
    meta = fixture(name)["meta"]
    assert meta["ivcrtest_version"] == R_REFERENCE_VERSION
    assert meta["kind"] == kind


def test_grid_map_is_linear_in_delta(fixture):
    # Lemma 3.1: A has rows (r, -sqrt(1 - r^2)), and the grid map is A delta.
    for case in fixture("algebra")["cases"]:
        r, A = arr(case["r"]), arr(case["A"])
        np.testing.assert_allclose(A, np.column_stack([r, -np.sqrt(1 - r**2)]), atol=1e-15)
        np.testing.assert_allclose(cr_delta(case["rho"]), case["delta"], atol=1e-15)
        lam = A @ arr(case["delta"])
        np.testing.assert_allclose(lam, case["lambda"], atol=1e-15)
        np.testing.assert_allclose(lam, case["g_fun"], atol=1e-12)


def test_omega_is_the_delta_method_sandwich(fixture):
    for case in fixture("algebra")["cases"]:
        J, S = arr(case["jacobian"]), arr(case["Sigma_rho"])
        omega = J @ S @ J.T
        np.testing.assert_allclose((omega + omega.T) / 2, case["Omega"], rtol=1e-13)


def test_nearpd_output_is_positive_semidefinite(fixture):
    for case in fixture("nearpd")["cases"]:
        out = arr(case["output"])
        np.testing.assert_allclose(out, out.T, atol=1e-18)
        assert np.linalg.eigvalsh(out).min() > 0
        if case["name"] == "positive_definite":
            np.testing.assert_allclose(out, case["input"], atol=1e-17)


def test_residuals_match_r_lm(fixture):
    d = fixture("data")
    data = {k: arr(v) for k, v in d["data"].items()}
    for case in d["residualised"].values():
        rows = np.asarray(case["rows"]) - 1  # R indices are 1-based
        assert len(rows) == case["n"]
        X = np.column_stack([np.ones(len(rows))] + [data[h][rows] for h in d["controls"]])

        def resid(v):
            beta, *_ = np.linalg.lstsq(X, v, rcond=None)
            return v - X @ beta

        z = resid(data[case["instrument"]][rows])
        if case["flipped"]:
            z = -z
        np.testing.assert_allclose(resid(data["x"][rows]), case["x"], atol=1e-12)
        np.testing.assert_allclose(resid(data["y"][rows]), case["y"], atol=1e-12)
        np.testing.assert_allclose(z, case["z"], atol=1e-12)
        np.testing.assert_allclose(corr3(case["x"], case["y"], case["z"]), case["rho"], atol=1e-14)


def test_bootstrap_indices_reproduce_r(fixture):
    d = fixture("data")["residualised"]["z1"]
    x, y, z = arr(d["x"]), arr(d["y"]), arr(d["z"])
    cs = fixture("bootstrap")["cov_small"]
    idx = np.asarray(cs["indices"]) - 1
    assert idx.shape == (cs["B"], len(x))
    rb = np.array([corr3(x[i], y[i], z[i]) for i in idx])
    np.testing.assert_allclose(rb, cs["Rb"], atol=1e-14)
    s = np.cov(rb, rowvar=False)
    np.testing.assert_allclose((s + s.T) / 2, cs["S_raw"], rtol=1e-12)


def test_im_bootstrap_indices_reproduce_r(fixture):
    d = fixture("data")["residualised"]["z1"]
    x, y, z = arr(d["x"]), arr(d["y"]), arr(d["z"])
    im = fixture("bootstrap")["im_bootstrap"]
    r = arr(im["grid"])
    A = np.column_stack([r, -np.sqrt(1 - r**2)])
    idx = np.asarray(im["indices"]) - 1
    assert idx.shape == (im["B"], len(x))
    ends = []
    for i in idx:
        lam = A @ cr_delta(corr3(x[i], y[i], z[i]))
        ends.append([lam.min(), lam.max()])
    ends = np.array(ends)
    np.testing.assert_allclose(ends, im["result"]["boot_endpoints"], atol=1e-14)
    sd = ends.std(axis=0, ddof=1)
    np.testing.assert_allclose(sd, [im["result"]["se_l"], im["result"]["se_u"]], rtol=1e-12)


@pytest.mark.parametrize("route", ["reparam", "primitive"])
def test_cihybrid_draws_cover_every_bisection_step(fixture, route):
    case = next(c for c in fixture("cihybrid")["cases"] if c["route"] == route)
    s = case["settings"]
    k = len(case["deltahat"])
    draws = [arr(m) for m in case["draws"]]
    assert draws[0].shape == (s["Blarge"], k)
    assert all(m.shape == (s["B"], k) for m in draws[1:])
    # The bisection halves [0, z_{1 - alpha/2}] until it is narrower than tol.
    steps = math.ceil(math.log2(NormalDist().inv_cdf(1 - s["alpha"] / 2) / s["tol"]))
    assert len(draws) - 1 == steps
    # The large draw has the covariance it should, up to sampling error.
    np.testing.assert_allclose(np.cov(draws[0], rowvar=False).diagonal(),
                               arr(case["deltaSigma"]).diagonal(), rtol=0.25)


def test_union_bound_inputs_are_consistent(fixture):
    for case in fixture("union_bounds")["cases"]:
        al, sig = arr(case["Al"]), arr(case["deltaSigma"])
        v = al @ sig @ al.T
        np.testing.assert_allclose(np.sqrt(v.diagonal()), case["sigma_l"], rtol=1e-12)
        s = np.sqrt(v.diagonal())
        np.testing.assert_allclose(v / np.outer(s, s), case["corr_l"], atol=1e-12)
        for p in case["ciproj_p"]:
            assert p["p"] is None or 0 <= p["p"] <= 1 + case["eta"]


def test_card_mcub_intervals_match_the_paper(fixture):
    card = fixture("end_to_end")["card"]
    if card is None:
        pytest.skip("the fixtures were made without the wooldridge package")
    # Table 2 of Dzhumashev and Tursunalieva (2026), as reproduced in the README.
    table2 = {"fatheduc": (-0.013, 0.274), "motheduc": (0.005, 0.273),
              "nearc2": (0.001, 0.075), "nearc4": (-0.012, 0.095)}
    for z, expected in table2.items():
        np.testing.assert_allclose(card["instruments"][z]["mcub"], expected, atol=5e-4 + 1e-12)


def test_simulated_decisions(fixture):
    res = fixture("end_to_end")["sim"]["results"]
    # z1 and z3 are valid instruments in the simulation; z2 is correlated with u.
    assert res["zero_in_mcub"] == [True, False, True]
