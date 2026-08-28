"""Two-beam Takagi-Taupin for one frozen phase profile -> I(eta).

Everything is dimensionless:

- depth      ``zeta = z / L`` in [0, 1]
- coupling   ``alpha = L / Lambda_ext``  (Pendelloesung periods across the slab)
- rocking    ``eta``, the deviation parameter already multiplied by L

The common forward susceptibility ``chi_0`` is gauged out exactly: in symmetric
geometry it enters both equations with the same coefficient, so a shared phase
substitution removes it. Complex ``chi_0`` only leaves a common ``exp(-mu z)``,
which cancels in every normalized observable.

``alpha`` is real, which is the non-absorbing assumption: the coupling is really
``C sqrt(chi_h chi_hbar)`` and that reduces to ``C |chi_h|`` only when
``chi_hbar = chi_h*``. See DESIGN.md for the full assumption list.
``ponytail: real alpha drops the Borrmann effect. If mode-selective damping turns
out to matter in P8, make alpha complex -- the RK4 integrators already carry
complex state, so only the signature and the input validation change.``

Laue (transmission). Both beams travel +z; entrance ``D0(0)=1, Dh(0)=0``::

    dD0/dzeta = i pi alpha e^{+i phi} Dh
    dDh/dzeta = i pi alpha e^{-i phi} D0 + i eta Dh

Bragg (reflection). ``Dh`` travels -z, so ``d/ds_h = -d/dzeta``, and the boundary
conditions sit at opposite faces (``D0(0)=1``, ``Dh(1)=0``) -- a two-point BVP.
Rather than shoot, integrate the Riccati variable ``R = Dh / D0`` backward from
``R(1) = 0``, which turns it into an IVP in the stable direction::

    dR/dzeta = -i pi alpha e^{-i phi} - i eta R - i pi alpha e^{+i phi} R^2

Analytic ground truth, asserted in ``tests/test_core.py``:

- Laue, ``phi=0``, ``eta=0``  ->  ``|Dh(1)|^2 = sin^2(pi alpha)``
- Laue, real ``eta``          ->  ``|D0|^2 + |Dh|^2 = 1`` (lossless)
- Bragg, ``phi=0``            ->  total-reflection plateau for ``|eta| < 2 pi alpha``
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = ["PhaseField", "ZERO_PHASE", "required_nz", "solve_laue", "solve_bragg", "rocking"]


class PhaseField(Protocol):
    """Frozen phase profile ``phi(zeta)``, vectorized over depth."""

    def __call__(self, zeta: NDArray[np.float64]) -> NDArray[np.float64]: ...


def ZERO_PHASE(zeta):
    """Unmodulated reference."""
    return np.zeros_like(np.asarray(zeta, dtype=np.float64))


def required_nz(alpha: float, alpha_ac: float, per_cycle: int = 20) -> int:
    """Depth steps needed to resolve both oscillations in the integrand.

    The phase grating turns over ``alpha_ac`` times across the slab and the
    Pendelloesung solution turns over ``alpha`` times; RK4 needs a handful of
    steps per cycle of the faster one.
    """
    return max(256, int(per_cycle * max(alpha, alpha_ac)))


def _phase_factor(phi, zeta: float) -> NDArray[np.complex128]:
    return np.exp(1j * np.asarray(phi(np.float64(zeta)), dtype=np.float64))


def solve_laue(
    alpha: float,
    eta: ArrayLike = 0.0,
    phi: PhaseField | None = None,
    n_z: int = 512,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Integrate the Laue system from ``zeta=0`` to ``1`` (fixed-step RK4).

    Vectorized over ``eta``: the phase profile is shared, so a whole rocking
    curve costs one sweep. Returns the exit amplitudes ``(D0, Dh)``.
    """
    if alpha < 0:
        raise ValueError(f"alpha must be non-negative, got {alpha}")
    if n_z < 1:
        raise ValueError(f"n_z must be >= 1, got {n_z}")

    eta_a = np.asarray(eta, dtype=np.float64)
    D0 = np.ones(eta_a.shape, dtype=np.complex128)
    Dh = np.zeros(eta_a.shape, dtype=np.complex128)
    phase = ZERO_PHASE if phi is None else phi
    c = 1j * np.pi * alpha
    h = 1.0 / n_z

    def f(zeta, d0, dh):
        e = _phase_factor(phase, zeta)
        return c * e * dh, c * np.conj(e) * d0 + 1j * eta_a * dh

    for i in range(n_z):
        z = i * h
        k1a, k1b = f(z, D0, Dh)
        k2a, k2b = f(z + 0.5 * h, D0 + 0.5 * h * k1a, Dh + 0.5 * h * k1b)
        k3a, k3b = f(z + 0.5 * h, D0 + 0.5 * h * k2a, Dh + 0.5 * h * k2b)
        k4a, k4b = f(z + h, D0 + h * k3a, Dh + h * k3b)
        D0 = D0 + (h / 6.0) * (k1a + 2 * k2a + 2 * k3a + k4a)
        Dh = Dh + (h / 6.0) * (k1b + 2 * k2b + 2 * k3b + k4b)

    return D0, Dh


def solve_bragg(
    alpha: float,
    eta: ArrayLike = 0.0,
    phi: PhaseField | None = None,
    n_z: int = 512,
) -> NDArray[np.complex128]:
    """Integrate the Bragg Riccati equation backward from ``R(1)=0`` to ``zeta=0``.

    Returns the complex amplitude reflectivity ``R(0) = Dh(0)/D0(0)`` at the
    entrance surface; the measured reflectivity is ``|R|^2``.
    """
    if alpha < 0:
        raise ValueError(f"alpha must be non-negative, got {alpha}")
    if n_z < 1:
        raise ValueError(f"n_z must be >= 1, got {n_z}")

    eta_a = np.asarray(eta, dtype=np.float64)
    R = np.zeros(eta_a.shape, dtype=np.complex128)
    phase = ZERO_PHASE if phi is None else phi
    c = 1j * np.pi * alpha
    h = -1.0 / n_z  # backward

    def f(zeta, r):
        e = _phase_factor(phase, zeta)
        return -c * np.conj(e) - 1j * eta_a * r - c * e * r * r

    for i in range(n_z):
        z = 1.0 + i * h
        k1 = f(z, R)
        k2 = f(z + 0.5 * h, R + 0.5 * h * k1)
        k3 = f(z + 0.5 * h, R + 0.5 * h * k2)
        k4 = f(z + h, R + h * k3)
        R = R + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    if not np.all(np.isfinite(R)):
        raise FloatingPointError("Bragg Riccati integration diverged; raise n_z or lower alpha")
    return R


def rocking(
    alpha: float,
    eta: ArrayLike,
    phi: PhaseField | None = None,
    geometry: str = "laue",
    n_z: int = 512,
) -> NDArray[np.float64]:
    """Diffracted intensity ``I(eta)`` for one frozen snapshot.

    ``geometry='laue'`` returns ``|Dh(L)|^2``; ``'bragg'`` returns ``|R(0)|^2``.
    """
    if geometry == "laue":
        _, Dh = solve_laue(alpha, eta, phi=phi, n_z=n_z)
        return np.abs(Dh) ** 2
    if geometry == "bragg":
        return np.abs(solve_bragg(alpha, eta, phi=phi, n_z=n_z)) ** 2
    raise ValueError(f"geometry must be 'laue' or 'bragg', got {geometry!r}")
