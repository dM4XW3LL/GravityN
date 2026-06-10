"""
body.py
───────
Python-side data model for one body in the N-body simulator.
 
This module is GUI-aware (it holds references to matplotlib artists) but
has no dependency on tkinter, so it can also be used by headless scripts
that only need to track body metadata alongside the physics engine.
 
The actual physics state (x, y, vx, vy, …) lives inside the C
``CSimulation`` struct managed by :class:`~nbody_engine.NBodyEngine`.
A ``Body`` instance does NOT duplicate that state; it stores only the
parameters used to *spawn* the body (needed for JSON save/load) and the
matplotlib artist handles used to *draw* it.

See difference between Body class and CBody class in documentation.
"""

from __future__ import annotations

# ─────────────────────────────────────────────────────────────────────────────
#  Internal counter — gives each body a unique ID within a session
# ─────────────────────────────────────────────────────────────────────────────
 
_BODY_COUNTER: int = 0

# ─────────────────────────────────────────────────────────────────────────────
#  Body
# ─────────────────────────────────────────────────────────────────────────────

class Body:
    """
    GUI-side representation of one simulated body.
 
    Attributes are split into three groups:
 
    **Identity / appearance**
        ``name``, ``mass``, ``color``, ``marker_size``
 
    **Spawn parameters** (stored so the system can be saved to JSON and
    reloaded later).  One group is populated depending on ``spawn_mode``;
    the other group's attributes are ``None``.
 
    - *Keplerian* (``spawn_mode == "elements"``) —
      ``a``, ``e``, ``M0_deg``, ``M_star``
    - *Cartesian* (``spawn_mode == "cartesian"``) —
      ``x0``, ``y0``, ``vx0``, ``vy0``
 
    **matplotlib artists** (set by the GUI when the body is added to the plot)
        ``dot``, ``trail_line``, ``label_text``, ``trail_x``, ``trail_y``
    """

    # ── Solar system presets ──────────────────────────────────────────────
    # Masses in solar masses; orbital elements from JPL mean elements (J2000).
    # M0_deg = 0 places every planet at periapsis at t = 0; this is a
    # simplification that produces visually correct orbits without requiring
    # real epoch data.
    PRESETS: dict[str, dict] = {
        "Mercury": dict(a=0.387,  e=0.206, M0_deg=0.0, mass=1.651e-7, color="#b5b5b5"),
        "Venus":   dict(a=0.723,  e=0.007, M0_deg=0.0, mass=2.447e-6, color="#e8cda0"),
        "Earth":   dict(a=1.000,  e=0.017, M0_deg=0.0, mass=3.003e-6, color="#4fa3e0"),
        "Mars":    dict(a=1.524,  e=0.093, M0_deg=0.0, mass=3.213e-7, color="#c1440e"),
        "Jupiter": dict(a=5.203,  e=0.049, M0_deg=0.0, mass=9.545e-4, color="#c88b3a"),
        "Saturn":  dict(a=9.537,  e=0.057, M0_deg=0.0, mass=2.858e-4, color="#e4d191"),
        "Uranus":  dict(a=19.19,  e=0.047, M0_deg=0.0, mass=4.365e-5, color="#7de8e8"),
        "Neptune": dict(a=30.07,  e=0.009, M0_deg=0.0, mass=5.149e-5, color="#5b7fde"),
    }

    def __init__(
        self,
        name: str,
        mass: float,
        color: str,
        spawn_mode: str = "elements",
        marker_size: int = 7,
        # Keplerian spawn parameters
        a:       float | None = None,
        e:       float | None = None,
        M0_deg:  float | None = None,
        M_star:  float | None = None,
        # Cartesian spawn parameters
        x:  float | None = None,
        y:  float | None = None,
        vx: float | None = None,
        vy: float | None = None,
    ) -> None:
        global _BODY_COUNTER
        _BODY_COUNTER += 1

        # ── Identity ──────────────────────────────────────────────────────
        self.id = _BODY_COUNTER
        self.name        = name
        self.mass        = mass
        self.color       = color
        self.marker_size = marker_size

        # ── Spawn mode ────────────────────────────────────────────────────
        #: Either ``"elements"`` (Keplerian) or ``"cartesian"``.
        self.spawn_mode  = spawn_mode

        # ── Keplerian spawn parameters ────────────────────────────────────
        self.a      = a        #: Semi-major axis (AU), or None.
        self.e      = e        #: Eccentricity, or None.
        self.M0_deg = M0_deg   #: Initial mean anomaly (degrees), or None.
        self.M_star = M_star   #: Central body mass (M☉), or None.

        # ── Cartesian spawn parameters ────────────────────────────────────
        self.x0  = x    #: Initial x position (AU), or None.
        self.y0  = y    #: Initial y position (AU), or None.
        self.vx0 = vx   #: Initial x velocity (AU yr⁻¹), or None.
        self.vy0 = vy   #: Initial y velocity (AU yr⁻¹), or None.

        # ── matplotlib artists (set by the GUI) ───────────────────────────
        #: The moving planet marker (``matplotlib.lines.Line2D``).
        self.dot        = None
        #: The fading position trail (``matplotlib.lines.Line2D``).
        self.trail_line = None
        #: The name label drawn near the dot (``matplotlib.text.Text``).
        self.label_text = None
 
        #: Trail buffer — x positions accumulated each animation tick.
        self.trail_x: list[float] = []
        #: Trail buffer — y positions accumulated each animation tick.
        self.trail_y: list[float] = []
    
    # ── Serialisation helpers ─────────────────────────────────────────────

    def to_dict(self) -> dict:
        """
        Serialise this body to a JSON-compatible dictionary.
 
        The dictionary contains enough information to reconstruct the body
        via :meth:`from_dict`.  matplotlib artists are not included.
        """
        entry: dict = {
            "name":       self.name,
            "mass":       self.mass,
            "color":      self.color,
            "spawn_mode": self.spawn_mode,
        }
        if self.spawn_mode == "elements":
            entry.update(a=self.a, e=self.e,
                         M0_deg=self.M0_deg, M_star=self.M_star)
        else:
            entry.update(x=self.x0, y=self.y0, vx=self.vx0, vy=self.vy0)
        return entry
    
    @classmethod
    def from_dict(cls, data: dict) -> "Body":
        """
        Reconstruct a :class:`Body` from a dictionary produced by
        :meth:`to_dict`.
 
        Raises:
            KeyError:   if a required field is missing.
            ValueError: if a numeric field cannot be parsed.
        """
        mode  = data.get("spawn_mode", "elements")
        name  = str(data["name"])
        mass  = float(data["mass"])
        color = str(data["color"])

        if mode == "elements":
            return cls(
                name=name, mass=mass, color=color, spawn_mode="elements",
                a=float(data["a"]), e=float(data["e"]),
                M0_deg=float(data["M0_deg"]), M_star=float(data["M_star"]),
            )
        else:
            return cls(
                name=name, mass=mass, color=color, spawn_mode="cartesian",
                x=float(data["x"]), y=float(data["y"]),
                vx=float(data["vx"]), vy=float(data["vy"]),
            )

    def __repr__(self) -> str:
        return f"Body({self.name!r}, mass={self.mass:.3e} M☉, mode={self.spawn_mode!r})"





