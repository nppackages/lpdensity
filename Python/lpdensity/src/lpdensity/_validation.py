import numpy as np
import math


def as_numeric_vector(value, name):
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise Exception(f"{name} should be numeric.") from None

    if array.ndim == 0:
        array = array.reshape(1)
    elif array.ndim == 1:
        array = array.reshape(-1)
    elif array.ndim == 2 and 1 in array.shape:
        array = array.reshape(-1)
    else:
        raise Exception(f"{name} should be one-dimensional.")

    if len(array) == 0:
        raise Exception(f"{name} cannot be empty.")
    if np.isnan(array).any():
        raise Exception(f"{name} should not contain missing values.")
    return array


def as_scalar(value, name):
    array = as_numeric_vector(value, name)
    if len(array) != 1:
        raise Exception(f"{name} should be a scalar.")
    return float(array[0])


def as_nonnegative_integer(value, name):
    try:
        scalar = as_scalar(value, name)
    except Exception:
        raise Exception(f"Option {name} incorrectly specified.") from None
    if math.ceil(scalar) < 0:
        raise Exception(f"Option {name} incorrectly specified.")
    return math.ceil(scalar)


def normalize_bwselect(value):
    aliases = {
        "mse": "mse-dpi",
        "imse": "imse-dpi",
        "rot": "mse-rot",
        "irot": "imse-rot",
    }
    if value is None:
        return "mse-dpi"
    if isinstance(value, str):
        text = value
    else:
        values = np.asarray(value).reshape(-1)
        text = "" if len(values) == 0 else str(values[0])
    text = text.lower()
    if len(text) == 0:
        return "mse-dpi"
    return aliases.get(text, text)
