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
FONT_TITLE = ("Courier New", 14, "bold")

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

# ─────────────────────────────────────────────────────────────────────────────
#  Main application window
# ─────────────────────────────────────────────────────────────────────────────

class NBodySimApp(tk.Tk):
    """
    Main application window for the N-body gravitational simulator.
    Inspired by the previous Keplerian Orbital Simulator layout.
    Three-column layout:
 
    +--------------+-------------------+-------------+
    | Left panel   | Matplotlib canvas | Info panel  |
    | (controls)   | (live simulation) | (collapsible|
    +--------------+-------------------+-------------+
 
    The left panel contains body-spawning controls, simulation controls
    (play/pause, integrator, speed), and the active-body list.
 
    The info panel shows live physical quantities for the selected body
    and system-wide conservation diagnostics.
    """


    # Days of simulation time advanced per animation tick.
    DAYS_PER_TICK = 2

    # How many days of position trail to keep behind each body.
    TRAIL_DAYS = 365*3

    # Tkinter ``after()`` interval between animation frames (ms).
    ANIM_INTERVAL: int = 30

    def __init__(self, lib) -> None:
        super().__init__()
        self.engine         = NBodyEngine(lib)
        self._user_data_dir = _get_user_data_dir()

        self.configure(bg=BG_DARK)
        self.resizable(True, True)
        self.title("GravityN - N-body Gravitational Simulator")

        # ── Simulation state ──────────────────────────────────────────────
        self.bodies: list[Body] = []
        self.running = False
        self._anim_id = None
        self._integrator = "leapfrog"

        # Reference energy recorded when the first body is added; used to
        # compute the relative energy drift shown in the info panel.
        self._E0: float | None = None

        # ── GUI state ─────────────────────────────────────────────────────
        self._next_color = ACCENT
        self._selected_body: Body | None = None
        self._selected_idx: int | None = None
        self._info_visible = False
        self._spawn_mode = tk.StringVar(value="elements")

        self._build_ui()
        self._style_ttk()
        self._draw_sun()

        self.after(100, self._resize_check)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ─────────────────────────────────────────────────────────────────────────────
    #  UI construction
    # ─────────────────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(0, weight=0)
        self.rowconfigure(0, weight=1)

        self._build_left_panel()
        self._build_canvas()
        self._build_info_panel()

    # ── Left panel ────────────────────────────────────────────────────────

    def _build_left_panel(self) -> None:
        panel = tk.Frame(self,bg=BG_PANEL, width=310)
        panel.grid(row=0,column=0,sticky="nsew")
        panel.grid_propagate(False)
        panel.columnconfigure(0, weight=1)

        row = 0

        #Title
        tk.Label(
            panel, text="N-BODY SIM",
            font=FONT_TITLE,
            fg=ACCENT, bg=BG_PANEL, pady=10
        ).grid(row=row, column=0, sticky="ew"); row+=1
        tk.Frame(panel,bg=BORDER,height=1).grid(
            row=row, column=0, sticky="ew"); row +=1
        
        # ── Preset selector ───────────────────────────────────────────────
        sec = self._section(panel, "PRESETS", row); row+=1
        self._preset_var = tk.StringVar(value="── Custom ──")
        preset_names = (["── Custom ──", "⊙ Full Solar System"]
                        + list(Body.PRESETS.keys()))
        
        self._preset_combo = ttk.Combobox(
            sec, textvariable=self._preset_var,
            values=preset_names, state="readonly", font=FONT_MONO
        )
        self._preset_combo.pack(fill="x", padx=8, pady=(0, 6))
        self._preset_combo.bind("<<ComboboxSelected>>",
                                self._on_preset_selected)
        
        # ── Spawn mode toggle ─────────────────────────────────────────────
        sec2 = self._section(panel, "SPAWN MODE", row); row += 1
        mode_row = tk.Frame(sec2, bg=BG_PANEL)
        mode_row.pack(fill="x", padx=8, pady=4)
 
        for text, value in [("Orbital elements", "elements"),
                             ("Cartesian  (x, y, vx, vy)", "cartesian")]:
            tk.Radiobutton(
                mode_row, text=text, variable=self._spawn_mode,
                value=value, command=self._on_spawn_mode_changed,
                bg=BG_PANEL, fg=TEXT_MAIN, selectcolor=BG_WIDGET,
                activebackground=BG_PANEL, font=FONT_LABEL
            ).pack(side="left", padx=(0, 8))

        # ── Orbital elements section ──────────────────────────────────────
        self._sec_elements = self._section(panel, "ORBITAL PARAMETERS", row)
        row += 1
        self._params: dict[str, tk.StringVar] = {}
 
        for label, key, default in [
            ("Name",                         "name",    "MyBody"),
            ("Mass (M\u2609)",               "mass",    "3.003e-6"),
            ("Semi-major axis (AU)",         "a",       "1.0"),
            ("Eccentricity",                 "e",       "0.0"),
            ("M\u2080 \u2014 init. mean anomaly (\u00b0)", "M0_deg", "0.0"),
            ("Central mass (M\u2609)",       "M_star",  "1.0"),
        ]:
            self._add_param_row(self._sec_elements, label, key, default)
 
        tk.Label(
            self._sec_elements,
            text=" \u2191  period auto-fills: P = a^1.5 / M_star^0.5",
            fg=TEXT_DIM, bg=BG_PANEL, font=FONT_SMALL,
            wraplength=270, justify="left"
        ).pack(anchor="w", padx=8, pady=(2, 6))

        # ── Cartesian section ─────────────────────────────────────────────
        self._sec_cartesian = self._section(panel, "CARTESIAN STATE", row)
        row += 1
 
        for label, key, default in [
            ("Name",           "c_name", "MyBody"),
            ("Mass (M\u2609)", "c_mass", "3.003e-6"),
            ("x  (AU)",        "c_x",    "1.0"),
            ("y  (AU)",        "c_y",    "0.0"),
            ("vx  (AU/yr)",    "c_vx",   "0.0"),
            ("vy  (AU/yr)",    "c_vy",   "6.283"),
        ]:
            self._add_param_row(self._sec_cartesian, label, key, default)
 
        # Cartesian section hidden by default
        self._sec_cartesian.grid_remove()

        # ── Color picker ──────────────────────────────────────────────────
        color_frame = tk.Frame(panel, bg=BG_PANEL)
        color_frame.grid(row=row, column=0, sticky="ew", padx=8, pady=2)
        row += 1
        tk.Label(color_frame, text="Color", fg=TEXT_DIM, bg=BG_PANEL,
                 font=FONT_LABEL, width=20, anchor="w").pack(side="left")
        self._color_btn = tk.Button(
            color_frame, bg=self._next_color, width=4,
            relief="flat", cursor="hand2", command=self._pick_color
        )
        self._color_btn.pack(side="left", padx=4)
 
        # ── Spawn button ──────────────────────────────────────────────────
        tk.Button(
            panel, text="  \uff0b  SPAWN BODY",
            font=("Courier New", 10, "bold"),
            bg=GREEN, fg=BG_DARK, activebackground="#2ea043",
            relief="flat", cursor="hand2", pady=6,
            command=self._spawn_body
        ).grid(row=row, column=0, sticky="ew", padx=10, pady=(8, 4))
        row += 1
 
        # ── Simulation controls ───────────────────────────────────────────
        sec3 = self._section(panel, "SIMULATION", row); row += 1
 
        ctrl = tk.Frame(sec3, bg=BG_PANEL)
        ctrl.pack(fill="x", padx=8, pady=4)
 
        self._play_btn = tk.Button(
            ctrl, text="\u25b6  PLAY", font=FONT_MONO,
            bg=ACCENT, fg=BG_DARK, activebackground="#79c0ff",
            relief="flat", cursor="hand2", width=10,
            command=self._toggle_animation
        )
        self._play_btn.pack(side="left", padx=(0, 4))
 
        tk.Button(
            ctrl, text="\u27f3  RESET", font=FONT_MONO,
            bg=BG_WIDGET, fg=TEXT_MAIN, activebackground=BORDER,
            relief="flat", cursor="hand2", width=10,
            command=self._reset_simulation
        ).pack(side="left")

        # Integrator selector
        integ_row = tk.Frame(sec3, bg=BG_PANEL)
        integ_row.pack(fill="x", padx=8, pady=(2, 4))
        tk.Label(integ_row, text="Integrator", fg=TEXT_DIM, bg=BG_PANEL,
                 font=FONT_LABEL).pack(side="left")
        self._integ_var = tk.StringVar(value="Leapfrog")
        ttk.Combobox(
            integ_row, textvariable=self._integ_var,
            values=["Leapfrog", "RK4"], state="readonly",
            font=FONT_SMALL, width=10
        ).pack(side="left", padx=6)
        self._integ_var.trace_add("write", self._on_integrator_changed)
 
        # Speed slider
        spd_row = tk.Frame(sec3, bg=BG_PANEL)
        spd_row.pack(fill="x", padx=8, pady=(2, 2))
        tk.Label(spd_row, text="Speed (days/tick)", fg=TEXT_DIM,
                 bg=BG_PANEL, font=FONT_LABEL).pack(side="left")
        self._speed_var = tk.IntVar(value=self.DAYS_PER_TICK)
        tk.Scale(
            spd_row, from_=1, to=100, orient="horizontal",
            variable=self._speed_var,
            bg=BG_PANEL, fg=TEXT_MAIN, troughcolor=BG_WIDGET,
            highlightthickness=0, bd=0, font=FONT_SMALL,
            showvalue=True, length=120
        ).pack(side="left", padx=6)
 
        # Trail slider
        trail_row = tk.Frame(sec3, bg=BG_PANEL)
        trail_row.pack(fill="x", padx=8, pady=(0, 6))
        tk.Label(trail_row, text="Trail (days)", fg=TEXT_DIM,
                 bg=BG_PANEL, font=FONT_LABEL).pack(side="left")
        self._trail_var = tk.IntVar(value=self.TRAIL_DAYS)
        tk.Scale(
            trail_row, from_=0, to=365 * 10, orient="horizontal",
            variable=self._trail_var,
            bg=BG_PANEL, fg=TEXT_MAIN, troughcolor=BG_WIDGET,
            highlightthickness=0, bd=0, font=FONT_SMALL,
            showvalue=True, length=120
        ).pack(side="left", padx=6)
 
        # Time display
        self._time_label = tk.Label(
            sec3, text="t = 0.0000 yr",
            fg=ACCENT, bg=BG_PANEL, font=("Courier New", 11, "bold")
        )
        self._time_label.pack(pady=(0, 6))

        # ── Active bodies list ────────────────────────────────────────────
        sec4 = self._section(panel, "ACTIVE BODIES", row); row += 1
 
        list_frame = tk.Frame(sec4, bg=BG_PANEL)
        list_frame.pack(fill="both", expand=True, padx=8, pady=(0, 6))
        scrollbar = tk.Scrollbar(list_frame, bg=BG_WIDGET)
        scrollbar.pack(side="right", fill="y")
 
        self._body_listbox = tk.Listbox(
            list_frame, font=FONT_SMALL,
            bg=BG_WIDGET, fg=TEXT_MAIN,
            selectbackground=ACCENT, selectforeground=BG_DARK,
            relief="flat", bd=0, highlightthickness=0,
            yscrollcommand=scrollbar.set, height=7,
        )
        self._body_listbox.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=self._body_listbox.yview)
        self._body_listbox.bind("<<ListboxSelect>>", self._on_body_selected)
 
        self._info_toggle_btn = tk.Button(
            sec4, text="\u2139  SHOW BODY INFO", font=FONT_SMALL,
            bg=BG_WIDGET, fg=ACCENT, activebackground=BORDER,
            relief="flat", cursor="hand2",
            command=self._toggle_info_panel
        )
        self._info_toggle_btn.pack(fill="x", padx=8, pady=(0, 4))
 
        tk.Button(
            sec4, text="\u2715  REMOVE SELECTED", font=FONT_SMALL,
            bg=BG_WIDGET, fg=ACCENT2, activebackground=BORDER,
            relief="flat", cursor="hand2",
            command=self._remove_selected_body
        ).pack(fill="x", padx=8, pady=(0, 4))

        # JSON save / load
        json_row = tk.Frame(sec4, bg=BG_PANEL)
        json_row.pack(fill="x", padx=8, pady=(0, 8))
        tk.Button(
            json_row, text="\U0001f4be  SAVE", font=FONT_SMALL,
            bg=BG_WIDGET, fg=GREEN, activebackground=BORDER,
            relief="flat", cursor="hand2", width=10,
            command=self._save_presets
        ).pack(side="left", padx=(0, 4))
        tk.Button(
            json_row, text="\U0001f4c2  LOAD", font=FONT_SMALL,
            bg=BG_WIDGET, fg=ACCENT, activebackground=BORDER,
            relief="flat", cursor="hand2", width=10,
            command=self._load_presets
        ).pack(side="left")
 
        # Footer
        tk.Frame(panel, bg=BORDER, height=1).grid(
            row=row, column=0, sticky="ew"); row += 1
        tk.Label(
            panel,
            text="Engine: C (nbody.c + kepler.c)  |  GUI: Python/Tk",
            fg=TEXT_DIM, bg=BG_PANEL, font=("Courier New", 7)
        ).grid(row=row, column=0, pady=6)



 





