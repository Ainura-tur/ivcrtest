"""Correlation-restriction (CR) test for instrument validity.

This is a Python port, in progress, of the R package ivcrtest
(https://github.com/Ainura-tur/ivcrtest). It implements the test of

    Dzhumashev, R. and A. Tursunalieva (2026). A test for instrumental
    variable validity using a partially identified correlation restriction.
    Econometric Reviews. https://doi.org/10.1080/07474938.2026.2741548

The port is checked against reference results recorded from the R package;
see python/README.md. No functions are available yet.
"""

__version__ = "0.1.0.dev0"

#: Version of the R package whose results this port is tested against.
R_REFERENCE_VERSION = "0.2.2"

__all__ = ["__version__", "R_REFERENCE_VERSION"]
