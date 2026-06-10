"""
nbody_engine.py
───────────────
Python interface to the compiled N-body C library (nbody.so / .dll / .dylib).

This module is intentionally GUI-free.  It can be imported from a tkinter
application, a Jupyter notebook, a command-line script, or a test or anything
that needs to drive the N-body physics engine from Python.

Contents:
  - load_nbody_lib()   — locate and load the shared library
  - CBody              — ctypes mirror of the C `Body` struct
  - CSimulation        — ctypes mirror of the C `Simulation` struct
  - NBodyEngine        — ergonomic wrapper around the C API

Physical units (identical to nbody.c / kepler.c):
  Length — Astronomical Unit  (AU)
  Time   — Julian year        (yr)
  Mass   — Solar mass         (M☉)
  G      = 4π²  AU³ yr⁻² M☉⁻¹
"""

import ctypes
import math
import os
import platform
import sys


# ─────────────────────────────────────────────────────────────────────────────
#  Constants — must match the #define values in nbody.h exactly
# ─────────────────────────────────────────────────────────────────────────────

#: Maximum number of bodies in one simulation (mirrors NBODY_MAX_BODIES).
NBODY_MAX_BODIES: int = 64

#: Maximum length of a body name string including the NUL terminator
#: (mirrors NBODY_NAME_LEN).
NBODY_NAME_LEN: int = 64


# ─────────────────────────────────────────────────────────────────────────────
#  1. Library loader
# ─────────────────────────────────────────────────────────────────────────────

def load_nbody_lib() -> ctypes.CDLL:
    """
    Locate and load the compiled nbody shared library.

    Search order:
      1. The directory containing this script  (src/gui/).
      2. The project root  (two levels up from src/gui/).
      3. The current working directory.

    This allows the app to be launched from any of:
      - ``python3 src/gui/nbodysim.py``   (project root)
      - ``python3 nbodysim.py``           (from inside src/gui/)

    Returns:
        A :class:`ctypes.CDLL` instance with all function signatures
        pre-registered.

    Raises:
        FileNotFoundError: if no matching library file is found in any of
        the search directories.
    """
    # Check if running inside a PyInstaller executable bundle
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        search_dirs = [sys._MEIPASS]
    else:
        script_dir   = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
        search_dirs  = [script_dir, project_root, os.getcwd()]

    system = platform.system()
    if system == "Windows":
        candidates = ["nbody.dll"]
    elif system == "Darwin":
        candidates = ["nbody.dylib", "nbody.so"]
    else:
        candidates = ["nbody.so", "nbody.dylib"]

    for directory in search_dirs:
        for name in candidates:
            path = os.path.join(directory, name)
            if os.path.exists(path):
                lib = ctypes.CDLL(path)
                _setup_lib_signatures(lib)
                return lib

    raise FileNotFoundError(
        f"Could not find the N-body shared library.\n"
        f"Expected one of: {candidates}\n"
        f"Searched in: {search_dirs}\n\n"
        f"Compile with (from the project root):\n"
        f"  Linux/macOS: gcc -O2 -fPIC -shared -o nbody.so "
        f"src/core/nbody.c src/core/kepler.c -lm\n"
        f"  Windows:     gcc -O2 -fPIC -shared -o nbody.dll "
        f"src/core/nbody.c src/core/kepler.c"
    )


def _setup_lib_signatures(lib: ctypes.CDLL) -> None:
    """
    Register the C type signatures for every function called via ctypes.

    This is mandatory: without explicit argtypes / restype declarations,
    ctypes assumes all arguments are C ints and returns are C ints, which
    silently corrupts double values on most platforms.

    Args:
        lib: The loaded CDLL object to annotate in place.
    """
    dbl  = ctypes.c_double
    pdbl = ctypes.POINTER(ctypes.c_double)

    # ── Lifecycle ─────────────────────────────────────────────────────────
    lib.nbody_init.argtypes  = [ctypes.c_void_p]
    lib.nbody_init.restype   = None

    lib.nbody_free.argtypes  = [ctypes.c_void_p]
    lib.nbody_free.restype   = None

    lib.nbody_add_body.argtypes  = [ctypes.c_void_p, ctypes.c_void_p]
    lib.nbody_add_body.restype   = ctypes.c_int

    lib.nbody_remove_body.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.nbody_remove_body.restype  = ctypes.c_int

    # ── Initial conditions ────────────────────────────────────────────────
    lib.elements_to_state.argtypes = [dbl, dbl, dbl, dbl,
                                       pdbl, pdbl, pdbl, pdbl]
    lib.elements_to_state.restype  = None

    # ── Core physics ──────────────────────────────────────────────────────
    lib.nbody_compute_accelerations.argtypes = [ctypes.c_void_p]
    lib.nbody_compute_accelerations.restype  = None

    # ── Integrators ───────────────────────────────────────────────────────
    lib.nbody_step_leapfrog.argtypes = [ctypes.c_void_p, dbl]
    lib.nbody_step_leapfrog.restype  = None

    lib.nbody_step_rk4.argtypes = [ctypes.c_void_p, dbl]
    lib.nbody_step_rk4.restype  = None

    # ── Diagnostics ───────────────────────────────────────────────────────
    lib.nbody_total_energy.argtypes        = [ctypes.c_void_p]
    lib.nbody_total_energy.restype         = dbl

    lib.nbody_total_momentum.argtypes      = [ctypes.c_void_p, pdbl, pdbl]
    lib.nbody_total_momentum.restype       = None

    lib.nbody_barycentre.argtypes          = [ctypes.c_void_p, pdbl, pdbl]
    lib.nbody_barycentre.restype           = None

    lib.nbody_body_kinetic_energy.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.nbody_body_kinetic_energy.restype  = dbl


# ─────────────────────────────────────────────────────────────────────────────
#  2. ctypes struct mirrors  (must match nbody.h layout exactly)
# ─────────────────────────────────────────────────────────────────────────────

class CBody(ctypes.Structure):
    """
    Mirror of the C ``Body`` struct from nbody.h.

    Field order and types must match the C declaration byte-for-byte.
    ctypes handles alignment automatically for native types.
    Verified size: 120 bytes on all supported platforms.

    .. code-block:: c

        typedef struct {
            double x, y, vx, vy, ax, ay;
            double mass;
            char   name[64];
        } Body;
    """
    _fields_ = [
        ("x",    ctypes.c_double),
        ("y",    ctypes.c_double),
        ("vx",   ctypes.c_double),
        ("vy",   ctypes.c_double),
        ("ax",   ctypes.c_double),
        ("ay",   ctypes.c_double),
        ("mass", ctypes.c_double),
        ("name", ctypes.c_char * NBODY_NAME_LEN),
    ]


class CSimulation(ctypes.Structure):
    """
    Mirror of the C ``Simulation`` struct from nbody.h.

    An instance of this struct *is* the complete simulation state.
    Python owns the memory; the C engine reads and writes into it via
    pointer.  Verified size: 7696 bytes on all supported platforms.

    .. code-block:: c

        typedef struct {
            Body   bodies[64];
            int    n;
            double t;
        } Simulation;
    """
    _fields_ = [
        ("bodies", CBody * NBODY_MAX_BODIES),
        ("n",      ctypes.c_int),
        ("t",      ctypes.c_double),
    ]


# ─────────────────────────────────────────────────────────────────────────────
#  NBodyEngine
# ─────────────────────────────────────────────────────────────────────────────

class NBodyEngine:
    """
    Ergonomic Python interface to the N-body C library.

    Owns a single :class:`CSimulation` struct and exposes methods that map
    to the C API.  All ctypes boilerplate and unit conversions live
    here so the rest of the application deals only with plain Python floats.

    Example::

        lib    = load_nbody_lib()
        engine = NBodyEngine(lib)

        engine.add_body_elements("Earth", mass=3.003e-6,
                                 a=1.0, e=0.017, M0_deg=0.0, M_star=1.0)
        for _ in range(1000):
            engine.step(1.0 / 1000)

        x, y = engine.position(0)
    """

    #: Unit conversion: 1 AU yr⁻¹ = 4.74047 km s⁻¹  (same as kepler.c)
    _AU_YR_TO_KM_S: float = 4.74047

    def __init__(self, lib: ctypes.CDLL) -> None:
        self._lib = lib
        self._sim = CSimulation()
        lib.nbody_init(ctypes.byref(self._sim))

    # ── Lifecycle ─────────────────────────────────────────────────────────

    def reset(self) -> None:
        """Remove all bodies and reset simulation time to zero."""
        self._lib.nbody_free(ctypes.byref(self._sim))
        self._lib.nbody_init(ctypes.byref(self._sim))

    # ── Body management ───────────────────────────────────────────────────

    def add_body_cartesian(self, name: str, mass: float,
                           x: float, y: float,
                           vx: float, vy: float) -> bool:
        """
        Add a body with a Cartesian initial state.

        After adding, accelerations are recomputed so that the very first
        leapfrog half-kick uses the correct forces.

        Args:
            name: Human-readable label (truncated to 63 characters).
            mass: Gravitational mass in M☉.  Must be > 0.
            x, y:   Initial position in AU.
            vx, vy: Initial velocity in AU yr⁻¹.

        Returns:
            ``True`` on success, ``False`` if the simulation is already
            at capacity (:data:`NBODY_MAX_BODIES` bodies).
        """
        b      = CBody()
        b.x    = x;    b.y  = y
        b.vx   = vx;   b.vy = vy
        b.mass = mass
        b.name = name.encode("utf-8")[:NBODY_NAME_LEN - 1]

        result = self._lib.nbody_add_body(ctypes.byref(self._sim),
                                          ctypes.byref(b))
        if result == 0:
            self._lib.nbody_compute_accelerations(ctypes.byref(self._sim))
            return True
        return False

    def add_body_elements(self, name: str, mass: float,
                          a: float, e: float,
                          M0_deg: float, M_star: float) -> bool:
        """
        Add a body from Keplerian orbital elements.

        Calls the C ``elements_to_state`` function to convert the elements
        to a Cartesian state vector, then delegates to
        :meth:`add_body_cartesian`.

        Args:
            name:    Human-readable label.
            mass:    Gravitational mass in M☉.
            a:       Semi-major axis in AU.  Must be > 0.
            e:       Eccentricity.  Must satisfy 0 ≤ e < 1.
            M0_deg:  Initial mean anomaly in **degrees**.
            M_star:  Mass of the central body in M☉.  Must be > 0.

        Returns:
            ``True`` on success, ``False`` if the simulation is full.
        """
        x  = ctypes.c_double(0.0)
        y  = ctypes.c_double(0.0)
        vx = ctypes.c_double(0.0)
        vy = ctypes.c_double(0.0)

        self._lib.elements_to_state(
            ctypes.c_double(a),
            ctypes.c_double(e),
            ctypes.c_double(math.radians(M0_deg)),
            ctypes.c_double(M_star),
            ctypes.byref(x), ctypes.byref(y),
            ctypes.byref(vx), ctypes.byref(vy),
        )
        return self.add_body_cartesian(name, mass,
                                       x.value, y.value,
                                       vx.value, vy.value)

    def remove_body(self, idx: int) -> bool:
        """
        Remove the body at zero-based index ``idx``.

        Remaining bodies shift down to fill the gap, preserving their
        relative order (matching the behaviour of ``nbody_remove_body`` in C).
        Accelerations are recomputed after removal.

        Returns:
            ``True`` on success, ``False`` if ``idx`` is out of range.
        """
        result = self._lib.nbody_remove_body(ctypes.byref(self._sim),
                                              ctypes.c_int(idx))
        if result == 0:
            self._lib.nbody_compute_accelerations(ctypes.byref(self._sim))
            return True
        return False

    # ── Integration ───────────────────────────────────────────────────────

    def step(self, dt_years: float, integrator: str = "leapfrog") -> None:
        """
        Advance the simulation by one timestep.

        Args:
            dt_years:   Timestep in years.  Must be > 0.
            integrator: ``"leapfrog"`` (default, symplectic) or ``"rk4"``
                        (4th-order, non-symplectic).  Any other string falls
                        back to leapfrog.
        """
        if integrator == "rk4":
            self._lib.nbody_step_rk4(ctypes.byref(self._sim),
                                      ctypes.c_double(dt_years))
        else:
            self._lib.nbody_step_leapfrog(ctypes.byref(self._sim),
                                           ctypes.c_double(dt_years))

    # ── State readback ────────────────────────────────────────────────────

    @property
    def n(self) -> int:
        """Number of active bodies currently in the simulation."""
        return self._sim.n

    @property
    def t(self) -> float:
        """Current simulation time in years."""
        return self._sim.t

    def position(self, idx: int) -> tuple[float, float]:
        """Return ``(x, y)`` position in AU for body at index ``idx``."""
        b = self._sim.bodies[idx]
        return b.x, b.y

    def velocity(self, idx: int) -> tuple[float, float, float]:
        """
        Return ``(vx, vy, speed)`` for body at index ``idx``.

        ``vx`` and ``vy`` are in AU yr⁻¹; ``speed`` is the scalar
        magnitude converted to km s⁻¹.
        """
        b   = self._sim.bodies[idx]
        spd = math.sqrt(b.vx**2 + b.vy**2) * self._AU_YR_TO_KM_S
        return b.vx, b.vy, spd

    def mass(self, idx: int) -> float:
        """Return the mass of body ``idx`` in M☉."""
        return self._sim.bodies[idx].mass

    def name(self, idx: int) -> str:
        """Return the name of body ``idx`` as a Python string."""
        return self._sim.bodies[idx].name.decode("utf-8")

    def distance_between(self, i: int, j: int) -> float:
        """Return the distance between bodies ``i`` and ``j`` in AU."""
        bi, bj = self._sim.bodies[i], self._sim.bodies[j]
        dx = bj.x - bi.x
        dy = bj.y - bi.y
        return math.sqrt(dx * dx + dy * dy)

    # ── Diagnostics ───────────────────────────────────────────────────────

    def total_energy(self) -> float:
        """
        Total mechanical energy of the system in M☉ AU² yr⁻².

        For a leapfrog-integrated isolated system this should oscillate
        around a constant.  A secular drift indicates too large a timestep.
        """
        return self._lib.nbody_total_energy(ctypes.byref(self._sim))

    def total_momentum(self) -> tuple[float, float]:
        """
        Total linear momentum ``(px, py)`` in M☉ AU yr⁻¹.

        Should be conserved to machine precision by any integrator.
        """
        px = ctypes.c_double(0.0)
        py = ctypes.c_double(0.0)
        self._lib.nbody_total_momentum(ctypes.byref(self._sim),
                                        ctypes.byref(px), ctypes.byref(py))
        return px.value, py.value

    def barycentre(self) -> tuple[float, float]:
        """
        Position of the system barycentre (centre of mass) in AU.

        Returns ``(cx, cy)``.
        """
        cx = ctypes.c_double(0.0)
        cy = ctypes.c_double(0.0)
        self._lib.nbody_barycentre(ctypes.byref(self._sim),
                                    ctypes.byref(cx), ctypes.byref(cy))
        return cx.value, cy.value

    def body_kinetic_energy(self, idx: int) -> float:
        """
        Kinetic energy of body ``idx`` in M☉ AU² yr⁻².

        Equal to ½ m v².  Does not include gravitational potential energy.
        """
        return self._lib.nbody_body_kinetic_energy(ctypes.byref(self._sim),
                                                    ctypes.c_int(idx))