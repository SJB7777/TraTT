"""DynamicFTT -- transient formation dynamics of acoustic Bragg satellites.

The solver produces one thing, the time-resolved rocking curve ``I(eta, t)``.
Every result in the paper is an analysis of that cube; no observable is added to
the solver itself.

Layout::

    core.py     frozen phi(zeta) -> I(eta)      (Laue IVP, Bragg Riccati)
    drive.py    phi(zeta; t)                    (traveling wave + causal envelope)
    cube.py     one parameter point -> I(eta,t)
    observe.py  I(eta,t) -> p_n(t), eta_c(t), sigma(t) -> formation times
    sweep.py    (alpha, S) grid -> t_f/t_tr map

See DESIGN.md for the physics conventions and the analytic ground truth.
"""

from dftt.core import required_nz, rocking, solve_bragg, solve_laue
from dftt.drive import Drive

__all__ = ["Drive", "required_nz", "rocking", "solve_bragg", "solve_laue"]
