#!/usr/bin/env python3
"""
nbodysim.py
────────────
Main entry point for the N-body Gravitational Simulator GUI.
 
This file is responsible purely for the user interface: building the tkinter
window, driving the matplotlib animation loop, and connecting user actions
(spawning bodies, changing integrators, toggling the info panel) to the
physics engine.
 
All C interop lives in :mod:`nbody_engine`.
All body data modelling lives in :mod:`body`.
 
Usage (from the project root)::
 
    gcc -O2 -fPIC -shared -o nbody.so src/core/nbody.c src/core/kepler.c -lm
    python3 src/gui/nbodysim.py
"""

import sys
import os
import platform
import ctypes
import math
import tkinter as tk
from tkinter import ttk, messagebox, colorchooser, filedialog
import json
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ─────────────────────────────────────────────────────────────────────────────
#  1.  Load the C shared library
# ─────────────────────────────────────────────────────────────────────────────

def load_nbody_lib():
    """
    Locate and load the compiled nbody shared library.
 
    Search order: the directory containing this script, then the project
    root (two levels up), then the current working directory.  This mirrors
    the behaviour of orbitsim.py so both apps can be launched the same way.
 
    Raises FileNotFoundError if no matching library is found.
    """

    #Check if running inside a PyInstaller executable bundle
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        search_dirs = [sys._MEIPASS]
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.abspath(os.path.join(script_dir,"..",".."))
        search_dirs = [script_dir, project_root, os.getcwd()]

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

def _setup_lib_signatures(lib):
    """
    Register the C type signatures for every function we call via ctypes.
 
    This is mandatory: without explicit argtypes / restype declarations,
    ctypes assumes all arguments are C ints and returns are C ints, which
    silently corrupts double values on most platforms.
    """
    dbl = ctypes.c_double
    pdbl = ctypes.POINTER(ctypes.c_double)
    pint = ctypes.POINTER(ctypes.c_int)

    # ── Simulation struct size helper (not a C function — Python only) ────
    # We pass the Simulation struct by pointer, so we need its size.
    # The struct layout is defined by nbody.h:
    #   Body  = 6 doubles (x,y,vx,vy,ax,ay) + 1 double (mass) + 64 chars (name)
    #         = 7*8 + 64 = 120 bytes
    #   Simulation = Body[64] + int (n) + double (t)
    #              = 64*120 + 8 + 8 = 7696 bytes  (+ possible padding)
    # We use ctypes structures so Python manages the memory correctly.

    # ── Lifecycle ─────────────────────────────────────────────────────────
    lib.nbody_init.argtypes = [ctypes.c_void_p]
    lib.nbody_init.restype = None

    lib.nbody_free.argtypes = [ctypes.c_void_p]
    lib.nbody_free.restype = None

    lib.nbody_add_body.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    lib.nbody_add_body.restype = ctypes.c_int

    lib.nbody_remove_body.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.nbody_remove_body.restype = ctypes.c_int

    # ── Initial conditions ────────────────────────────────────────────────
    lib.elements_to_state.argtypes = [dbl, dbl, dbl, dbl, pdbl, pdbl, pdbl, pdbl]
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
    lib.nbody_total_energy.argtypes       = [ctypes.c_void_p]
    lib.nbody_total_energy.restype        = dbl
 
    lib.nbody_total_momentum.argtypes     = [ctypes.c_void_p, pdbl, pdbl]
    lib.nbody_total_momentum.restype      = None
 
    lib.nbody_barycentre.argtypes         = [ctypes.c_void_p, pdbl, pdbl]
    lib.nbody_barycentre.restype          = None
 
    lib.nbody_body_kinetic_energy.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.nbody_body_kinetic_energy.restype  = dbl


# ─────────────────────────────────────────────────────────────────────────────
#  2.  ctypes struct mirrors  (must match nbody.h exactly)
# ─────────────────────────────────────────────────────────────────────────────

NBODY_MAX_BODIES = 64
NBODY_NAME_LEN   = 64


class CBody(ctypes.Structure):
    """
    Mirror of the C `Body` struct from nbody.h.
 
    Field order and types must match the C declaration byte-for-byte.
    ctypes handles alignment automatically for native types.
    """
    _fields_ = [
        ("x", ctypes.c_double),
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
    Mirror of the C `Simulation` struct from nbody.h.
 
    An instance of this struct IS the simulation state.  We pass a pointer
    to it into every C function.  Python owns the memory; C writes into it.
    """
    _fields_ = [
        ("bodies", CBody * NBODY_MAX_BODIES),
        ("n",      ctypes.c_int),
        ("t",      ctypes.c_double),
    ]

