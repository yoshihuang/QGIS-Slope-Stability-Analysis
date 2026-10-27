# -*- coding: utf-8 -*-
"""Material definitions and shear-strength models.

Supported strength models
-------------------------
``mc``          Mohr-Coulomb, effective stress:  tau = c' + sigma_n' tan(phi')
``undrained``   Undrained (phi = 0):             tau = Su   (Su given in ``cohesion``)
``hb``          Generalised Hoek-Brown (2002) rock mass criterion, converted to an
                instantaneous (tangent) c_i / phi_i at the current normal stress
                (Balmer 1952 relations).
``none``        No strength (e.g. water, very weak fill) - tau = 0.

Units: kPa, kN/m3 and degrees throughout (any consistent unit set works if the
unit weight of water is changed accordingly).
"""

import math
import numpy as np

STRENGTH_MODELS = ("mc", "undrained", "hb", "none")

_MODEL_ALIASES = {
    "": "mc", "mc": "mc", "mohr": "mc", "mohr-coulomb": "mc", "mohr_coulomb": "mc",
    "drained": "mc", "cphi": "mc",
    "undrained": "undrained", "su": "undrained", "cu": "undrained", "phi0": "undrained",
    "hb": "hb", "hoek": "hb", "hoek-brown": "hb", "hoek_brown": "hb",
    "none": "none", "water": "none", "zero": "none",
}


def normalise_model(name):
    key = str(name or "").strip().lower()
    return _MODEL_ALIASES.get(key, "mc")


class Material(object):
    """A soil / rock material."""

    def __init__(self, name="Material", gamma=19.0, cohesion=0.0, phi=0.0,
                 gamma_sat=None, model="mc", sigci=0.0, gsi=0.0, mi=0.0,
                 disturbance=0.0, ru=0.0, use_piezo=True,
                 c_std=0.0, phi_std=0.0, gamma_std=0.0):
        self.name = str(name)
        self.gamma = float(gamma)
        self.gamma_sat = float(gamma_sat) if gamma_sat not in (None, "") else float(gamma)
        self.cohesion = float(cohesion)
        self.phi = float(phi)
        self.model = normalise_model(model)
        self.sigci = float(sigci or 0.0)
        self.gsi = float(gsi or 0.0)
        self.mi = float(mi or 0.0)
        self.disturbance = float(disturbance or 0.0)
        self.ru = float(ru or 0.0)
        self.use_piezo = bool(use_piezo)
        self.c_std = float(c_std or 0.0)
        self.phi_std = float(phi_std or 0.0)
        self.gamma_std = float(gamma_std or 0.0)
        self._hb_table = None
        if self.model == "hb":
            self._build_hb_table()

    # ------------------------------------------------------------------ utils
    def copy(self, **changes):
        m = Material(self.name, self.gamma, self.cohesion, self.phi, self.gamma_sat,
                     self.model, self.sigci, self.gsi, self.mi, self.disturbance,
                     self.ru, self.use_piezo, self.c_std, self.phi_std, self.gamma_std)
        for k, v in changes.items():
            setattr(m, k, v)
        if m.model == "hb":
            m._build_hb_table()
        return m

    @property
    def is_nonlinear(self):
        return self.model == "hb"

    def validate(self):
        """Return a list of (level, message) problems for pre-flight checks."""
        from .i18n import tr
        msgs = []
        if self.gamma <= 0:
            msgs.append(("error", tr("Unit weight must be greater than zero.")))
        if self.gamma_sat < self.gamma * 0.5:
            msgs.append(("warning", tr("Saturated unit weight is much lower than the unit weight.")))
        if self.model == "mc":
            if self.cohesion < 0:
                msgs.append(("error", tr("Cohesion cannot be negative.")))
            if not (0.0 <= self.phi < 90.0):
                msgs.append(("error", tr("Friction angle must be between 0 and 90 degrees.")))
            if self.cohesion == 0 and self.phi == 0:
                msgs.append(("warning", tr("Material has no shear strength (c = 0 and phi = 0).")))
        elif self.model == "undrained":
            if self.cohesion <= 0:
                msgs.append(("error", tr("Undrained strength Su (cohesion field) must be greater than zero.")))
        elif self.model == "hb":
            if self.sigci <= 0 or self.mi <= 0 or not (0 < self.gsi <= 100):
                msgs.append(("error", tr("Hoek-Brown requires sigci > 0, mi > 0 and 0 < GSI <= 100.")))
            if not (0.0 <= self.disturbance <= 1.0):
                msgs.append(("error", tr("Hoek-Brown disturbance factor D must be between 0 and 1.")))
        if not (0.0 <= self.ru < 1.0):
            msgs.append(("error", tr("Pore pressure ratio ru must be between 0 and 1.")))
        return msgs

    # --------------------------------------------------------- Hoek-Brown
    def hb_parameters(self):
        gsi, mi, d = self.gsi, self.mi, self.disturbance
        mb = mi * math.exp((gsi - 100.0) / (28.0 - 14.0 * d))
        s = math.exp((gsi - 100.0) / (9.0 - 3.0 * d))
        a = 0.5 + (math.exp(-gsi / 15.0) - math.exp(-20.0 / 3.0)) / 6.0
        return mb, s, a

    def _build_hb_table(self):
        if self.sigci <= 0 or self.mi <= 0 or self.gsi <= 0:
            self._hb_table = None
            return
        mb, s, a = self.hb_parameters()
        sci = self.sigci
        s3_min = -s * sci / mb * 0.999          # uniaxial tensile limit
        # dense near the tensile limit, geometric towards high stress
        s3 = np.concatenate([
            np.linspace(s3_min, 0.0, 60, endpoint=False),
            np.geomspace(sci * 1e-6, sci * 2.0, 240),
        ])
        base = np.maximum(mb * s3 / sci + s, 1e-12)
        s1 = s3 + sci * base ** a
        ds1ds3 = 1.0 + a * mb * base ** (a - 1.0)
        sn = (s1 + s3) / 2.0 - (s1 - s3) / 2.0 * (ds1ds3 - 1.0) / (ds1ds3 + 1.0)
        tau = (s1 - s3) * np.sqrt(ds1ds3) / (ds1ds3 + 1.0)
        order = np.argsort(sn)
        sn, tau = sn[order], tau[order]
        keep = np.concatenate([[True], np.diff(sn) > 1e-12])
        sn, tau = sn[keep], tau[keep]
        slope = np.gradient(tau, sn)
        self._hb_table = (sn, tau, slope)

    # ------------------------------------------------------- strength
    def linear_parameters(self, sigma_n=None):
        """Return (c, tan_phi) - tangent parameters at normal stress ``sigma_n``.

        ``sigma_n`` is only used by non-linear models; it may be a scalar or array.
        """
        if self.model == "mc":
            return self.cohesion, math.tan(math.radians(self.phi))
        if self.model == "undrained":
            return self.cohesion, 0.0
        if self.model == "none":
            return 0.0, 0.0
        # Hoek-Brown
        if self._hb_table is None:
            return 0.0, 0.0
        sn_t, tau_t, slope_t = self._hb_table
        sn = 0.0 if sigma_n is None else sigma_n
        sn = np.clip(np.asarray(sn, dtype=float), sn_t[0], sn_t[-1])
        tau = np.interp(sn, sn_t, tau_t)
        tphi = np.maximum(np.interp(sn, sn_t, slope_t), 0.0)
        c = np.maximum(tau - sn * tphi, 0.0)
        if np.ndim(c) == 0:
            return float(c), float(tphi)
        return c, tphi
