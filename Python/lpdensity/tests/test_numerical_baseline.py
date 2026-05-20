from pathlib import Path

import numpy as np
import pandas as pd

from lpdensity import lpbwdensity, lpdensity
from lpdensity import lpbwdensity_fn

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _baseline_data():
    index = np.arange(200, dtype=float)
    return pd.Series(np.linspace(0.05, 3.95, 200) + 0.05 * np.sin(index), name="v1")


def _weights():
    index = np.arange(1, 201, dtype=float)
    return 1 + 0.1 * (index % 5), 1 + 0.05 * (index % 7)


def test_lpdensity_broad_numerical_baseline():
    data = _baseline_data()
    cweights, pweights = _weights()
    expected = pd.read_csv(FIXTURES / "lpdensity-baseline-broad.csv")
    cases = {
        "fixed_triangular_density": {"bw": 0.5, "grid": np.arange(0, 4.5, 0.5)},
        "kernel_uniform_density": {"bw": 0.5, "grid": np.arange(0, 5), "kernel": "uniform"},
        "kernel_epanechnikov_density": {"bw": 0.5, "grid": np.arange(0, 5), "kernel": "epanechnikov"},
        "cdf_triangular": {"bw": 0.75, "grid": np.arange(0, 5), "p": 1, "q": 2, "v": 0},
        "weighted_triangular_density": {
            "bw": 0.65,
            "grid": np.arange(0, 5),
            "Cweights": cweights,
            "Pweights": pweights,
        },
    }

    for case, kwargs in cases.items():
        observed = lpdensity(data, **kwargs).Estimate
        case_expected = expected.loc[expected["case"] == case, observed.columns]
        np.testing.assert_allclose(
            observed.to_numpy(dtype=float),
            case_expected.to_numpy(dtype=float),
            rtol=1e-10,
            atol=1e-10,
            err_msg=case,
        )


def test_lpbwdensity_broad_numerical_baseline():
    data = _baseline_data()
    expected = pd.read_csv(FIXTURES / "lpbwdensity-baseline-broad.csv")

    for selector in ["mse-dpi", "imse-dpi", "mse-rot", "imse-rot"]:
        case = f"bw_{selector}"
        observed = lpbwdensity(data, grid=np.arange(0, 5), bwselect=selector).table
        case_expected = expected.loc[expected["case"] == case, observed.columns]
        np.testing.assert_allclose(
            observed.to_numpy(dtype=float),
            case_expected.to_numpy(dtype=float),
            rtol=1e-6,
            atol=1e-6,
            err_msg=case,
        )


def test_lpdensity_unique_uses_last_occurrence_indices():
    observed = lpbwdensity_fn.__lpdensityUnique(np.array([1.0, 1.0, 2.0, 2.0, 2.0, 5.0]))
    np.testing.assert_array_equal(observed["unique"], np.array([1.0, 2.0, 5.0]))
    np.testing.assert_array_equal(observed["freq"], np.array([2, 3, 1]))
    np.testing.assert_array_equal(observed["index"], np.array([1, 4, 5]))

    tied = lpbwdensity_fn.__lpdensityUnique(np.array([4.0, 4.0, 4.0]))
    np.testing.assert_array_equal(tied["unique"], np.array([4.0]))
    np.testing.assert_array_equal(tied["freq"], np.array([3]))
    np.testing.assert_array_equal(tied["index"], np.array([2]))


def test_confint_uniform_invalid_simulation_count_falls_back_to_pointwise():
    fit = lpdensity(_baseline_data(), bw=[0.5], grid=np.arange(0, 4.5, 0.5))

    observed = fit.confint(CIuniform=True, CIsimul=1)
    expected = fit.confint(CIuniform=False)
    np.testing.assert_allclose(observed.to_numpy(), expected.to_numpy())

    np.random.seed(42)
    simulated = fit.confint(CIuniform=True, CIsimul=10)
    assert simulated.shape == (9, 7)
    assert np.isfinite(simulated.to_numpy()).all()


def test_documented_scalar_inputs_are_accepted():
    data = _baseline_data()

    scalar = lpdensity(data, bw=0.5, grid=0.5, scale=1).Estimate
    vector = lpdensity(data, bw=[0.5], grid=[0.5], scale=[1]).Estimate
    np.testing.assert_allclose(scalar.to_numpy(), vector.to_numpy(), rtol=1e-12, atol=1e-12)

    bandwidth = lpbwdensity(data, grid=0.5, bwselect="mse")
    assert bandwidth.bwselect == "mse-dpi"
    assert bandwidth.table.shape == (1, 4)
