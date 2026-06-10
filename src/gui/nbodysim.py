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

import json
import math
import os
import platform
import sys
import tkinter as tk
from tkinter import ttk, messagebox, colorchooser, filedialog
 
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from nbody_engine import load_nbody_lib, NBodyEngine
from body import Body


# ─────────────────────────────────────────────────────────────────────────────
#  1. Visual theme
# ─────────────────────────────────────────────────────────────────────────────

BG_DARK    = "#0d1117"
BG_PANEL   = "#161b22"
BG_WIDGET  = "#21262d"
BORDER     = "#30363d"
TEXT_MAIN  = "#e6edf3"
TEXT_DIM   = "#8b949e"
ACCENT     = "#58a6ff"
ACCENT2    = "#f78166"
SUN_COLOR  = "#ffe066"
GREEN      = "#3fb950"
RED        = "#f85149"
 
FONT_LABEL = ("Courier New", 9)
FONT_SMALL = ("Courier New", 8)
FONT_MONO  = ("Courier New", 10)


# ─────────────────────────────────────────────────────────────────────────────
#  Helper: user-data directory for JSON preset save / load
# ─────────────────────────────────────────────────────────────────────────────

def _get_user_data_dir() -> str:
    system = platform.system()
    if system == "Windows":
        base = os.environ.get("APPDATA", os.path.expanduser("~"))
        path = os.path.join(base, "NBodySim")
    else:
        path = os.path.join(os.path.expanduser("~"), ".nbody_sim")
    os.makedirs(path, exist_ok=True)
    return path




