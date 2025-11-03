"""
This module establishes helps the user input potential parameters for the
conical order they are expected in the rust code.

-CAM 2025
"""


def _get_or_zero(dict, key):
    v = dict.get(key)
    if v is None:
        return 0.0
    else:
        return v


def _make_canonical_form(dict):
    params = [
        "V",
        "r",
        "a",
        "W",
        "riv",
        "aiv",
        "WS",
        "ris",
        "ais",
        "SWS",
        "sri",
        "sai",
        "Vso",
        "rso",
        "aso",
        "rc",
    ]
    return [_get_or_zero(dict, p) for p in params]


def make_parameters(**kwargs):
    return _make_canonical_form(kwargs)
