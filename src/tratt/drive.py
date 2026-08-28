"""Traveling acoustic modulation -> frozen phase profile ``phi(zeta)``.

Slow time ``t = t_real / t_tr`` with wavefront transit ``t_tr = L / v_s``.
The drive of the abstract is::

    phi(zeta; t) = S * A(t - zeta) * sin(q*zeta - w*t)

with ``S = h.u0`` the modulation index, ``q = K*L``, ``w = Omega*t_tr``, and ``A``
a causal envelope of retarded time ``s = t - zeta`` (``A(s<0)=0`` reproduces the
wavefront cutoff ``z_max(t) = v_s t``).

**q and w are not independent.** For a non-dispersive traveling wave
``Omega = K*v_s``, so

    w = Omega*L/v_s = K*L = q = 2*pi*alpha_ac,     alpha_ac = L / Lambda_ac

leaving ``S`` and ``alpha_ac`` as the only drive parameters. With ``q = w`` the
phase field collapses to a pure function of retarded time,

    phi = -S * A(s) * sin(q*s),      s = t - zeta

i.e. for ``t >= 1`` the grating is a rigid profile sliding at ``v_s``, already in
its final form. That is why the formation time cannot exceed the transit time
under a hard step onset -- see DESIGN.md, fact 5.

The only knob that can push formation past the transit is a finite transducer
rise time ``rise = T_r / t_tr``, so it is an explicit parameter, not a smoothing
convenience.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = ["Drive"]


@dataclass(frozen=True)
class Drive:
    """Non-dispersive traveling acoustic modulation in dimensionless groups.

    Parameters
    ----------
    S:
        Modulation index ``h . u0``.
    alpha_ac:
        Acoustic thickness ``L / Lambda_ac``. Sets ``q = w = 2*pi*alpha_ac``, and
        hence the satellite spacing in ``eta``.
    rise:
        Transducer rise time in units of ``t_tr``. ``0`` = hard step onset
        (strictly causal); ``>0`` = C1 smoothstep over ``s in [0, rise]``.
    """

    S: float
    alpha_ac: float
    rise: float = 0.0

    def __post_init__(self) -> None:
        if self.alpha_ac <= 0:
            raise ValueError(f"alpha_ac must be positive, got {self.alpha_ac}")
        if self.rise < 0:
            raise ValueError(f"rise must be non-negative, got {self.rise}")

    @property
    def q(self) -> float:
        """Carrier spatial phase across the slab, ``K*L``. Equals ``w``."""
        return 2.0 * np.pi * self.alpha_ac

    @property
    def w(self) -> float:
        """Acoustic phase advanced per wavefront transit, ``Omega*t_tr``."""
        return self.q

    @property
    def period(self) -> float:
        """Acoustic period in slow-time units: ``2*pi/w = 1/alpha_ac``."""
        return 1.0 / self.alpha_ac

    def envelope(self, s: NDArray[np.float64]) -> NDArray[np.float64]:
        """Causal envelope of retarded time ``s`` (units of ``t_tr``)."""
        s = np.asarray(s, dtype=np.float64)
        if self.rise == 0.0:
            return np.where(s >= 0.0, 1.0, 0.0)
        x = np.clip(s / self.rise, 0.0, 1.0)
        return x * x * (3.0 - 2.0 * x)

    def phi(self, t: float):
        """Frozen phase profile ``phi(zeta)`` at slow time ``t``."""
        S, q = self.S, self.q

        def _phi(zeta):
            zeta = np.asarray(zeta, dtype=np.float64)
            return S * self.envelope(t - zeta) * np.sin(q * zeta - q * t)

        return _phi

    def eta_channel(self, n: int) -> float:
        """Center of diffraction channel ``n`` on the ``eta`` axis: ``eta_n = n*q``.

        In the kinematic limit ``Dh(1) = i pi alpha int_0^1 e^{-i phi(zeta)}
        e^{i eta (1-zeta)} dzeta``. Expanding
        ``e^{-i phi} = sum_n (-1)^n J_n(S) e^{i n (q zeta - q t)}`` leaves
        ``int_0^1 e^{i (n q - eta) zeta} dzeta``, which peaks at ``eta = n*q``.
        Positions are analytic; nothing is peak-fitted.

        ``n`` is the harmonic index of the traveling grating. Labeling the
        channels Stokes / anti-Stokes needs the frequency-resolved treatment --
        the stroboscopic model carries no photon-energy axis, so that labeling is
        deliberately not asserted here.

        For integer ``alpha_ac`` the finite-slab sinc of one channel vanishes
        exactly at every other channel center, which is why the Bessel ladder can
        be read off the peak heights with no leakage correction.
        """
        return n * self.q

    def times(self, t_max: float = 3.0, per_period: int = 32) -> NDArray[np.float64]:
        """Slow-time grid from 0 to ``t_max``, ``per_period`` samples per acoustic period."""
        if per_period < 1:
            raise ValueError("per_period must be >= 1")
        n = int(np.ceil(t_max / (self.period / per_period))) + 1
        return np.linspace(0.0, t_max, n)
