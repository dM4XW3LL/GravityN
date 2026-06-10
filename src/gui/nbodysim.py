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

    # ── Canvas ────────────────────────────────────────────────────────────
    def _build_canvas(self) -> None:
        frame = tk.Frame(self, bg=BG_DARK)
        frame.grid(row=0, column=1, sticky="nsew", padx=4, pady=4)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
 
        self.fig, self.ax = plt.subplots(figsize=(7, 7), facecolor=BG_DARK)
        self._setup_axes()
 
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self.canvas.draw()

    def _setup_axes(self) -> None:
        ax = self.ax
        ax.set_facecolor(BG_DARK)
        ax.tick_params(colors=TEXT_DIM, labelsize=7)
        for spine in ax.spines.values():
            spine.set_edgecolor(BORDER)
        ax.set_xlabel("x  (AU)", color=TEXT_DIM, fontsize=8,
                      fontfamily="monospace")
        ax.set_ylabel("y  (AU)", color=TEXT_DIM, fontsize=8,
                      fontfamily="monospace")
        ax.set_title("N-body Gravitational Simulator", color=TEXT_MAIN,
                     fontsize=10, fontfamily="monospace", pad=8)
        ax.set_aspect("equal")
        ax.grid(True, color=BORDER, linestyle="--", alpha=0.35,
                linewidth=0.5)
        ax.set_xlim(-12, 12)
        ax.set_ylim(-12, 12)

    def _draw_sun(self) -> None:
        """Draw the Sun marker at the origin with a layered glow effect."""
        self._sun_dot, = self.ax.plot(
            0, 0, "o", color=SUN_COLOR, markersize=12,
            zorder=10, markeredgewidth=0
        )
        for r, a in [(18, 0.06), (12, 0.10), (7, 0.18)]:
            self.ax.plot(0, 0, "o", color=SUN_COLOR, markersize=r,
                         alpha=a, zorder=9, markeredgewidth=0)
        self.canvas.draw_idle()


    # ── Info panel ────────────────────────────────────────────────────────
    def _build_info_panel(self) -> None:
        """Build the collapsible right-hand info panel."""
        self._info_frame = tk.Frame(self, bg=BG_PANEL, width=220)
        self._info_frame.columnconfigure(0, weight=1)
 
        tk.Label(
            self._info_frame, text="BODY  INFO",
            font=("Courier New", 11, "bold"),
            fg=ACCENT, bg=BG_PANEL, pady=8
        ).pack(fill="x")
        tk.Frame(self._info_frame, bg=BORDER, height=1).pack(fill="x")
 
        def info_row(parent: tk.Widget, label: str) -> tk.Label:
            """Return the value Label for one info row."""
            row = tk.Frame(parent, bg=BG_PANEL)
            row.pack(fill="x", padx=12, pady=2)
            tk.Label(row, text=label, fg=TEXT_DIM, bg=BG_PANEL,
                     font=FONT_LABEL, width=14, anchor="w").pack(side="left")
            val = tk.Label(row, text="\u2014", fg=TEXT_MAIN, bg=BG_PANEL,
                           font=FONT_MONO, anchor="w")
            val.pack(side="left")
            return val
 
        def section_label(parent: tk.Widget, text: str) -> None:
            tk.Label(parent, text=f" {text}", fg=TEXT_DIM, bg=BG_WIDGET,
                     font=("Courier New", 8, "bold"),
                     anchor="w", pady=6).pack(fill="x")
 
        def divider(parent: tk.Widget) -> None:
            tk.Frame(parent, bg=BORDER, height=1).pack(
                fill="x", padx=12, pady=4)
 
        # ── Body static fields ────────────────────────────────────────────
        section_label(self._info_frame, "BODY")
        self._info_name = info_row(self._info_frame, "Name")
        self._info_mass = info_row(self._info_frame, "Mass (M\u2609)")
        divider(self._info_frame)

        # ── Kinematic fields (updated every tick) ─────────────────────────
        section_label(self._info_frame, "KINEMATICS")
        self._info_x     = info_row(self._info_frame, "x")
        self._info_y     = info_row(self._info_frame, "y")
        self._info_speed = info_row(self._info_frame, "Speed")
        self._info_ke    = info_row(self._info_frame, "Kin. energy")
        divider(self._info_frame)

        # ── Distances (rebuilt dynamically each tick) ─────────────────────
        section_label(self._info_frame, "DISTANCES")
        self._info_dist_frame = tk.Frame(self._info_frame, bg=BG_PANEL)
        self._info_dist_frame.pack(fill="x", padx=12)
        divider(self._info_frame)
 
        # ── System-wide diagnostics ───────────────────────────────────────
        section_label(self._info_frame, "SYSTEM")
        self._info_E    = info_row(self._info_frame, "Total energy")
        self._info_dE   = info_row(self._info_frame, "Energy drift")
        self._info_p    = info_row(self._info_frame, "|momentum|")
        self._info_bary = info_row(self._info_frame, "Barycentre")
        divider(self._info_frame)

        # Placeholder shown when nothing is selected
        self._info_placeholder = tk.Label(
            self._info_frame,
            text="Click a body in the\nlist to inspect it.",
            fg=TEXT_DIM, bg=BG_PANEL,
            font=FONT_SMALL, justify="center", pady=10
        )
        self._info_placeholder.pack()

    # ─────────────────────────────────────────────────────────────────────────────
    #  Widget helpers
    # ─────────────────────────────────────────────────────────────────────────────

    def _section(self, parent: tk.Widget, title: str,
                 grid_row: int) -> tk.Frame:
        
        """Return a titled section frame, already grid-placed."""
        wrapper = tk.Frame(parent, bg=BG_PANEL)
        wrapper.grid(row=grid_row, column=0, sticky="ew")
        wrapper.columnconfigure(0, weight=1)
        hdr = tk.Frame(wrapper, bg=BG_WIDGET)
        hdr.pack(fill="x")
        tk.Label(hdr, text=f"  {title}",
                 font=("Courier New", 8, "bold"),
                 fg=TEXT_DIM, bg=BG_WIDGET,
                 anchor="w", pady=3).pack(fill="x")
        return wrapper

    def _add_param_row(self, parent: tk.Widget,
                       label: str, key: str, default: str) -> None:
        """Add one label + entry row to ``parent`` and register it in ``self._params``."""
        row = tk.Frame(parent, bg=BG_PANEL)
        row.pack(fill="x", padx=8, pady=2)
        tk.Label(row, text=label, fg=TEXT_DIM, bg=BG_PANEL,
                 font=FONT_LABEL, width=24, anchor="w").pack(side="left")
        var = tk.StringVar(value=default)
        self._params[key] = var
        tk.Entry(
            row, textvariable=var, font=FONT_MONO,
            bg=BG_WIDGET, fg=TEXT_MAIN, insertbackground=ACCENT,
            relief="flat", bd=4, width=10
        ).pack(side="left", fill="x", expand=True)

    def _style_ttk(self) -> None:
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("TCombobox",
            fieldbackground=BG_WIDGET, background=BG_WIDGET,
            foreground=TEXT_MAIN, selectbackground=ACCENT,
            selectforeground=BG_DARK, font=FONT_MONO,
            bordercolor=BORDER, arrowcolor=TEXT_DIM,
        )
    
    
    # ─────────────────────────────────────────────────────────────────────────────
    #  Spawn mode & Preset Handling
    # ─────────────────────────────────────────────────────────────────────────────

    def _on_spawn_mode_changed(self) -> None:
        if self._spawn_mode.get() == "elements":
            self._sec_cartesian.grid_remove()
            self._sec_elements.grid()
        else:
            self._sec_elements.grid_remove()
            self._sec_cartesian.grid()

    def _on_preset_selected(self, _event=None) -> None:
        name = self._preset_var.get()
        if name == "\u2500\u2500 Custom \u2500\u2500":
            return
        if name == "\u2299 Full Solar System":
            self._spawn_full_solar_system()
            return
        p = Body.PRESETS.get(name)
        if not p:
            return
        self._spawn_mode.set("elements")
        self._on_spawn_mode_changed()
        self._params["name"].set(name)
        self._params["mass"].set(str(p["mass"]))
        self._params["a"].set(str(p["a"]))
        self._params["e"].set(str(p["e"]))
        self._params["M0_deg"].set(str(p["M0_deg"]))
        self._params["M_star"].set("1.0")
        self._next_color = p["color"]
        self._color_btn.configure(bg=self._next_color)

    
    def _pick_color(self) -> None:
        color = colorchooser.askcolor(color=self._next_color,
                                       title="Choose body color")
        if color and color[1]:
            self._next_color = color[1]
            self._color_btn.configure(bg=self._next_color)
 
    def _on_integrator_changed(self, *_) -> None:
        self._integrator = ("rk4" if self._integ_var.get() == "RK4"
                            else "leapfrog")
        
    # ─────────────────────────────────────────────────────────────────────────────
    #  Body Spawning
    # ─────────────────────────────────────────────────────────────────────────────

    def _spawn_body(self) -> None:
        try:
            body = self._parse_and_build_body()
        except ValueError as exc:
            messagebox.showerror("Invalid Parameters", str(exc))
            return
        self._add_body_to_simulation(body)

    def _spawn_full_solar_system(self) -> None:
        """Clear existing bodies and spawn all 8 planets around the Sun."""
        if self.bodies:
            if not messagebox.askyesno(
                "Clear existing bodies?",
                "Spawning the full Solar System will clear all current "
                "bodies. Continue?"
            ):
                return
            self._reset_simulation()
 
        for planet_name, p in Body.PRESETS.items():
            body = Body(
                name=planet_name, mass=p["mass"], color=p["color"],
                spawn_mode="elements",
                a=p["a"], e=p["e"], M0_deg=p["M0_deg"], M_star=1.0,
            )
            self._add_body_to_simulation(body)
 
    def _parse_and_build_body(self) -> Body:
        """
        Read the active parameter fields and return a :class:`Body`.
 
        Raises:
            ValueError: with a descriptive message if any field is invalid.
        """
        errors: list[str] = []
 
        if self._spawn_mode.get() == "elements":
            name = self._params["name"].get().strip() or "Body"
            floats: dict[str, float] = {}
            for key, label in [("mass",   "Mass"),
                                ("a",      "Semi-major axis"),
                                ("e",      "Eccentricity"),
                                ("M0_deg", "M\u2080"),
                                ("M_star", "Central mass")]:
                try:
                    floats[key] = float(self._params[key].get())
                except ValueError:
                    errors.append(f"{label}: not a valid number")
 
            if "mass"   in floats and floats["mass"]   <= 0:
                errors.append("Mass must be > 0")
            if "a"      in floats and floats["a"]       <= 0:
                errors.append("Semi-major axis must be > 0")
            if "e"      in floats and not (0.0 <= floats["e"] < 1.0):
                errors.append("Eccentricity must be in [0, 1)")
            if "M_star" in floats and floats["M_star"] <= 0:
                errors.append("Central mass must be > 0")
 
            if errors:
                raise ValueError("\n".join(errors))
 
            return Body(
                name=name, mass=floats["mass"], color=self._next_color,
                spawn_mode="elements",
                a=floats["a"], e=floats["e"],
                M0_deg=floats["M0_deg"], M_star=floats["M_star"],
            )
 
        else:  # cartesian
            name = self._params["c_name"].get().strip() or "Body"
            floats = {}
            for key, label in [("c_mass", "Mass"), ("c_x", "x"),
                                ("c_y", "y"), ("c_vx", "vx"),
                                ("c_vy", "vy")]:
                try:
                    floats[key] = float(self._params[key].get())
                except ValueError:
                    errors.append(f"{label}: not a valid number")
 
            if "c_mass" in floats and floats["c_mass"] <= 0:
                errors.append("Mass must be > 0")
 
            if errors:
                raise ValueError("\n".join(errors))
 
            return Body(
                name=name, mass=floats["c_mass"], color=self._next_color,
                spawn_mode="cartesian",
                x=floats["c_x"], y=floats["c_y"],
                vx=floats["c_vx"], vy=floats["c_vy"],
            )
 
    def _add_body_to_simulation(self, body: Body) -> None:
        """
        Add a :class:`Body` to both the C engine and the matplotlib canvas.
 
        Reads back the initial position from the engine after adding so the
        trail buffer starts at the correct location.
        """
        if body.spawn_mode == "elements":
            ok = self.engine.add_body_elements(
                body.name, body.mass,
                body.a, body.e, body.M0_deg, body.M_star
            )
        else:
            ok = self.engine.add_body_cartesian(
                body.name, body.mass,
                body.x0, body.y0, body.vx0, body.vy0
            )
 
        if not ok:
            messagebox.showerror(
                "Simulation Full",
                f"Cannot add more than {NBodyEngine.__module__} bodies.")
            return
 
        idx = self.engine.n - 1
        x0, y0 = self.engine.position(idx)
 
        body.trail_line, = self.ax.plot(
            [], [], color=body.color, lw=1.2, alpha=0.55, zorder=3)
        body.dot, = self.ax.plot(
            [x0], [y0], "o", color=body.color,
            markersize=body.marker_size, zorder=6, markeredgewidth=0)
        body.label_text = self.ax.text(
            x0, y0 + 0.15, body.name,
            color=body.color, fontsize=6.5,
            fontfamily="monospace", zorder=7, alpha=0.85)
        body.trail_x = [x0]
        body.trail_y = [y0]
 
        self.bodies.append(body)
 
        # Set reference energy once we have at least two bodies
        # (one body alone has no potential energy to track).
        if self._E0 is None and self.engine.n > 1:
            self._E0 = self.engine.total_energy()
 
        self._refresh_body_list()
        self._update_axes_limits()
        self.canvas.draw_idle()

    # ─────────────────────────────────────────────────────────────────────────────
    #  Body Management
    # ─────────────────────────────────────────────────────────────────────────────

    def _remove_selected_body(self) -> None:
        sel = self._body_listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        self.engine.remove_body(idx)
 
        body = self.bodies[idx]
        for artist in (body.trail_line, body.dot, body.label_text):
            if artist is not None:
                artist.remove()
        self.bodies.pop(idx)
 
        if self._selected_idx == idx:
            self._selected_body = None
            self._selected_idx  = None
            self._info_placeholder.pack()
 
        # Recalibrate the reference energy after the body count changes.
        self._E0 = (self.engine.total_energy()
                    if self.engine.n > 1 else None)
 
        self._refresh_body_list()
        self._update_axes_limits()
        self.canvas.draw_idle()
 
    def _refresh_body_list(self) -> None:
        self._body_listbox.delete(0, tk.END)
        for b in self.bodies:
            self._body_listbox.insert(
                tk.END, f"  {b.name:<12}  {b.mass:.2e} M\u2609")
 
    def _update_axes_limits(self) -> None:
        """Auto-scale the plot to enclose all bodies with a 15% margin."""
        if not self.bodies:
            self.ax.set_xlim(-12, 12)
            self.ax.set_ylim(-12, 12)
            return
        reaches = []
        for b in self.bodies:
            if b.spawn_mode == "elements" and b.a is not None:
                reaches.append(b.a * (1 + b.e))
            elif b.x0 is not None:
                reaches.append(math.sqrt(b.x0**2 + b.y0**2) * 1.5)
        lim = max(max(reaches) * 1.15, 1.5) if reaches else 12.0
        self.ax.set_xlim(-lim, lim)
        self.ax.set_ylim(-lim, lim)
 
    def _on_body_selected(self, _event=None) -> None:
        sel = self._body_listbox.curselection()
        if not sel:
            self._selected_body = None
            self._selected_idx  = None
            self._info_placeholder.pack()
            return
        idx = sel[0]
        self._selected_body = self.bodies[idx]
        self._selected_idx  = idx
        self._info_placeholder.pack_forget()
        b = self._selected_body
        self._info_name.configure(text=b.name)
        self._info_mass.configure(text=f"{b.mass:.4e} M\u2609")

    # ─────────────────────────────────────────────────────────────────────────────
    #  Animation Loop
    # ─────────────────────────────────────────────────────────────────────────────
    
    def _toggle_animation(self) -> None:
        if self.running:
            self._stop_animation()
        else:
            self._start_animation()
 
    def _start_animation(self) -> None:
        self.running = True
        self._play_btn.configure(text="\u23f8  PAUSE", bg=ACCENT2)
        self._tick()
 
    def _stop_animation(self) -> None:
        self.running = False
        self._play_btn.configure(text="\u25b6  PLAY", bg=ACCENT)
        if self._anim_id:
            self.after_cancel(self._anim_id)
            self._anim_id = None
 
    def _tick(self) -> None:
        """
        One animation frame.
 
        Advances the C simulation by ``dt`` years, reads back all positions,
        updates matplotlib artists, and schedules the next frame.
        """
        if not self.running:
            return
 
        days_per_tick = self._speed_var.get()
        trail_days    = self._trail_var.get()
        dt_years      = days_per_tick / 365.25
        trail_years   = trail_days    / 365.25
 
        # Physics step
        self.engine.step(dt_years, self._integrator)
 
        # Update artists
        for i, body in enumerate(self.bodies):
            x, y = self.engine.position(i)
            body.trail_x.append(x)
            body.trail_y.append(y)
            body.dot.set_data([x], [y])
 
            # Trim trail to the requested time window
            while (len(body.trail_x) > 2 and
                   len(body.trail_x) * dt_years > trail_years):
                body.trail_x.pop(0)
                body.trail_y.pop(0)
 
            body.trail_line.set_data(body.trail_x, body.trail_y)
 
            if body.label_text is not None:
                body.label_text.set_position((x, y + 0.08))
 
        self._time_label.configure(text=f"t = {self.engine.t:.4f} yr")
        self._update_info_panel()
        self.canvas.draw_idle()
 
        self._anim_id = self.after(self.ANIM_INTERVAL, self._tick)
 
    def _update_info_panel(self) -> None:
        """
        Refresh all live fields in the info panel.
 
        Called every tick.  System diagnostics are always updated;
        per-body fields are only updated when a body is selected.
        """
        if not self._info_visible:
            return
 
        # System diagnostics
        if self.engine.n > 0:
            E      = self.engine.total_energy()
            px, py = self.engine.total_momentum()
            cx, cy = self.engine.barycentre()
            p_mag  = math.sqrt(px**2 + py**2)
 
            self._info_E.configure(text=f"{E:.4e}")
 
            if self._E0 is not None and self._E0 != 0.0:
                drift = abs((E - self._E0) / self._E0)
                color = (GREEN if drift < 1e-5
                         else ACCENT2 if drift < 1e-3
                         else RED)
                self._info_dE.configure(text=f"{drift:.2e}", fg=color)
            else:
                self._info_dE.configure(text="\u2014")
 
            self._info_p.configure(text=f"{p_mag:.4e}")
            self._info_bary.configure(
                text=f"({cx:.3f}, {cy:.3f}) AU")
 
        # Per-body fields
        if (self._selected_body is None
                or self._selected_idx is None
                or self._selected_idx >= self.engine.n):
            return
 
        idx         = self._selected_idx
        x, y        = self.engine.position(idx)
        _, _, spd   = self.engine.velocity(idx)
        ke          = self.engine.body_kinetic_energy(idx)
 
        self._info_x.configure(text=f"{x:.4f} AU")
        self._info_y.configure(text=f"{y:.4f} AU")
        self._info_speed.configure(text=f"{spd:.2f} km/s")
        self._info_ke.configure(
            text=f"{ke:.4e} M\u2609 AU\u00b2 yr\u207b\u00b2")
 
        # Distance rows — rebuilt each tick to handle bodies being added
        # or removed mid-simulation without stale entries.
        for widget in self._info_dist_frame.winfo_children():
            widget.destroy()
 
        for j, other in enumerate(self.bodies):
            if j == idx:
                continue
            dist = self.engine.distance_between(idx, j)
            row  = tk.Frame(self._info_dist_frame, bg=BG_PANEL)
            row.pack(fill="x", pady=1)
            tk.Label(row, text=f"\u2192 {other.name}", fg=TEXT_DIM,
                     bg=BG_PANEL, font=FONT_LABEL,
                     width=14, anchor="w").pack(side="left")
            tk.Label(row, text=f"{dist:.4f} AU", fg=TEXT_MAIN,
                     bg=BG_PANEL, font=FONT_MONO,
                     anchor="w").pack(side="left")


    # ─────────────────────────────────────────────────────────────────────────────
    #  Reset
    # ─────────────────────────────────────────────────────────────────────────────

    def _reset_simulation(self) -> None:
        """Remove all bodies, clear all artists, and reset the C engine."""
        was_running = self.running
        self._stop_animation()
 
        for body in self.bodies:
            for artist in (body.trail_line, body.dot, body.label_text):
                if artist is not None:
                    artist.remove()
 
        self.bodies.clear()
        self._selected_body = None
        self._selected_idx  = None
        self._E0            = None
 
        self.engine.reset()
        self._refresh_body_list()
 
        self.ax.cla()
        self._setup_axes()
        self._draw_sun()
        self._update_axes_limits()
        self.canvas.draw_idle()
        self._time_label.configure(text="t = 0.0000 yr")
 
        if was_running:
            self._start_animation()

    # ─────────────────────────────────────────────────────────────────────────────
    #  Info panel toggle
    # ─────────────────────────────────────────────────────────────────────────────

    def _toggle_info_panel(self) -> None:
        if self._info_visible:
            self._info_frame.grid_forget()
            self._info_visible = False
            self._info_toggle_btn.configure(
                text="\u2139  SHOW BODY INFO")
        else:
            self._info_frame.grid(row=0, column=2, sticky="nsew",
                                   padx=(0, 4), pady=4)
            self._info_frame.grid_propagate(False)
            self._info_visible = True
            self._info_toggle_btn.configure(
                text="\u2139  HIDE BODY INFO")
                
    # ─────────────────────────────────────────────────────────────────────────────
    #  JSON preset save / load
    # ─────────────────────────────────────────────────────────────────────────────

    def _save_presets(self) -> None:
        if not self.bodies:
            messagebox.showinfo("Nothing to save",
                                "Spawn at least one body first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save bodies", initialdir=self._user_data_dir,
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            initialfile="my_system.json",
        )
        if not path:
            return
        data = [b.to_dict() for b in self.bodies]
        with open(path, "w") as f:
            json.dump(data, f, indent=4)
        messagebox.showinfo("Saved",
                            f"Saved {len(data)} body/bodies to:\n{path}")
 
    def _load_presets(self) -> None:
        path = filedialog.askopenfilename(
            title="Load bodies", initialdir=self._user_data_dir,
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "r") as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as exc:
            messagebox.showerror("Load failed",
                                  f"Could not read file:\n{exc}")
            return
 
        spawned, errors = 0, []
        for i, entry in enumerate(data):
            try:
                body = Body.from_dict(entry)
                self._add_body_to_simulation(body)
                spawned += 1
            except (KeyError, ValueError) as exc:
                errors.append(
                    f"Entry {i+1} ({entry.get('name','?')}): {exc}")
 
        self._refresh_body_list()
        self._update_axes_limits()
        self.canvas.draw_idle()
 
        if errors:
            messagebox.showwarning(
                "Loaded with errors",
                f"Spawned {spawned} body/bodies.\n\nSkipped:\n"
                + "\n".join(errors))
        else:
            messagebox.showinfo("Loaded",
                                 f"Spawned {spawned} body/bodies.")

    # ─────────────────────────────────────────────────────────────────────────────
    #  Misc Functions (not yet categorized)
    # ─────────────────────────────────────────────────────────────────────────────

    def _resize_check(self) -> None:
        self.fig.set_facecolor(BG_DARK)
        self.after(2000, self._resize_check)
 
    def _on_close(self) -> None:
        self._stop_animation()
        self.destroy()
        sys.exit(0)


# ─────────────────────────────────────────────────────────────────────────────
#  Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    try:
        lib = load_nbody_lib()
    except FileNotFoundError as exc:
        root =tk.Tk()
        root.withdraw()
        messagebox.showerror("Library Not Found", str(exc))
        root.destroy()
        sys.exit(1)

    app = NBodySimApp(lib)
    app.geometry("1220x780")
    app.minsize(920,620)
    app.mainloop()

if __name__ == "__main__":
    main()
