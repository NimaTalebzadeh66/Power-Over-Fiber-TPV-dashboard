#!/usr/bin/env python3
"""Interactive thermal power-over-fiber dashboard.

Run: python fiber_thermal_gui.py
Install: pip install numpy pandas openpyxl matplotlib pillow

Keep this script, blackbody.xlsx, attenuation.xlsx, pof_scene.png and
pof_schematic.png together in one directory. No plots or data files are saved.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PV_BANDS = {"Si": 1100, "GaSb": 1700, "InGaAsSb": 2500}
FIBERS = {"Purchased": (0, 1), "UHV": (3, 4), "Low attenuation (OH)": (6, 7)}
TEMPERATURE_COLUMNS = {1500: 3, 2000: 2, 2500: 1}
EMISSIVE_W_CM2 = {1500: 28.7, 2000: 90.7, 2500: 221.5}
FIBER_COLORS = {"Purchased": "#ff918d", "UHV": "#a990ff", "Low attenuation (OH)": "#18d8e7"}
PV_COLORS = {"Si": "#f1be6f", "GaSb": "#a990ff", "InGaAsSb": "#70d9bd",
             "Custom": "#f596c8"}
BG = "#0e1230"
PANEL = "#171d40"
CARD = "#212a55"
TEXT = "#f0f2ff"
MUTED = "#9faed4"
CYAN = "#18d8e7"


def clean_xy(frame, x_column, y_column, name):
    x = pd.to_numeric(frame.iloc[1:, x_column], errors="coerce").to_numpy(float)
    y = pd.to_numeric(frame.iloc[1:, y_column], errors="coerce").to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 2 or np.any(x <= 0) or np.any(y < 0):
        raise ValueError(f"Invalid measurements for {name}")
    order = np.argsort(x)
    x, y = x[order], y[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError(f"Duplicate wavelengths in {name}")
    return x, y


def interpolate_with_endpoint_lines(x, y, grid):
    values = np.interp(grid, x, y)
    left, right = grid < x[0], grid > x[-1]
    values[left] = y[0] + (grid[left] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    values[right] = y[-1] + (grid[right] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return values


def gaussian_smooth(values, sigma_nm):
    if sigma_nm == 0:
        return values
    radius = int(np.ceil(4 * sigma_nm))
    offsets = np.arange(-radius, radius + 1)
    weights = np.exp(-0.5 * (offsets / sigma_nm) ** 2)
    weights /= weights.sum()
    return np.convolve(np.pad(values, radius, mode="edge"), weights, mode="valid")


def planck_exitance_nm(wavelength_nm, temperature):
    """Hemispherical blackbody spectral exitance, W cm^-2 nm^-1."""
    h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
    meters = np.asarray(wavelength_nm, dtype=float) * 1e-9
    exponent = h * c / (meters * k * temperature)
    return (2 * np.pi * h * c**2 / meters**5 / np.expm1(exponent)) * 1e-13


class SpectralModel:
    def __init__(self, directory=ROOT, smoothing_nm=20.0):
        self.wavelengths = np.arange(300.0, 2501.0, 1.0)
        blackbody = pd.read_excel(directory / "blackbody.xlsx", header=None)
        attenuation = pd.read_excel(directory / "attenuation.xlsx", header=None)
        self.sources = {}
        self.planck = {}
        self.alpha = {}
        self.measured_ranges = {}
        for temp, col in TEMPERATURE_COLUMNS.items():
            x, y = clean_xy(blackbody, 0, col, f"{temp} K")
            spectrum = interpolate_with_endpoint_lines(x, y, self.wavelengths)
            # The supplied samples match Planck's law. Use its shape only for
            # the unmeasured 2400–2500 nm tail, joined continuously at 2400 nm.
            tail = self.wavelengths > x[-1]
            if np.any(tail):
                spectrum[tail] = y[-1] * (
                    planck_exitance_nm(self.wavelengths[tail], temp)
                    / planck_exitance_nm(x[-1], temp)
                )
            self.sources[temp] = spectrum
            self.planck[temp] = planck_exitance_nm(self.wavelengths, temp)
            reference = planck_exitance_nm(x, temp)
            reference /= reference.max()
            difference = np.max(np.abs(reference - y))
            if difference > 0.02:
                raise ValueError(f"{temp} K source deviates from a Planck curve (max error {difference:.3f}); "
                                 "W/cm² calibration needs review")
        for name, cols in FIBERS.items():
            x, y = clean_xy(attenuation, *cols, name)
            self.measured_ranges[name] = (x[0], x[-1])
            values = np.maximum(interpolate_with_endpoint_lines(x, y, self.wavelengths), 0.0)
            self.alpha[name] = gaussian_smooth(values, smoothing_nm)

    def calculate(self, temperature, distance_km, bands, fibers):
        wl = self.wavelengths
        spectrum = self.sources[temperature]
        flux_spectrum = self.planck[temperature]
        total_planck = 5.670374419e-12 * temperature**4  # sigma*T^4 in W/cm²
        scaled_flux = EMISSIVE_W_CM2[temperature] / total_planck
        curves, results = {}, []
        for fiber in fibers:
            transmission = 10.0 ** (-self.alpha[fiber] * distance_km / 10.0)
            curves[fiber] = spectrum * transmission / np.max(spectrum)
            for label, low, high in bands:
                interior = wl[(wl > low) & (wl < high)]
                band_wl = np.concatenate(([low], interior, [high]))
                source_band = np.interp(band_wl, wl, spectrum)
                transmission_band = np.interp(band_wl, wl, transmission)
                flux_band = np.interp(band_wl, wl, flux_spectrum)
                incident = np.trapezoid(source_band, band_wl)
                remaining = np.trapezoid(source_band * transmission_band, band_wl)
                flux = np.trapezoid(flux_band * transmission_band, band_wl)
                results.append((label, fiber, remaining / incident, flux * scaled_flux))
        return spectrum / np.max(spectrum), curves, results


def start_gui(model):
    import tkinter as tk
    from tkinter import messagebox
    from PIL import Image, ImageOps, ImageTk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure

    class App:
        def __init__(self):
            self.root = tk.Tk()
            self.root.title("Power over Fiber | Thermal Link")
            self.root.geometry("1480x930")
            self.root.minsize(1130, 735)
            self.root.configure(bg=BG)
            self.temperature = tk.IntVar(value=2000)
            self.distance = tk.DoubleVar(value=1.0)
            self.pv_vars = {name: tk.BooleanVar(value=(name == "GaSb")) for name in PV_BANDS}
            self.pv_all_var = tk.BooleanVar(value=False)
            self.fiber_vars = {name: tk.BooleanVar(value=True) for name in FIBERS}
            self.custom_var = tk.BooleanVar(value=False)
            self.low_var = tk.StringVar(value="800")
            self.high_var = tk.StringVar(value="1800")
            self.pending = None
            self._build()
            self.schedule()

        def label(self, parent, text, size=11, color=TEXT, bold=False, **kw):
            return tk.Label(parent, text=text, bg=parent.cget("bg"), fg=color,
                            font=("Segoe UI", size, "bold" if bold else "normal"), **kw)

        def section(self, parent, title):
            shell = tk.Frame(parent, bg=BG)
            shell.pack(fill="x", pady=(0, 9))
            shell.pack_propagate(False)
            backdrop = tk.Canvas(shell, bg=BG, highlightthickness=0, bd=0)
            backdrop.place(x=0, y=0, relwidth=1, relheight=1)
            frame = tk.Frame(shell, bg=PANEL, padx=17, pady=12)
            window = backdrop.create_window(9, 4, anchor="nw", window=frame)

            def fit(event=None):
                width = max(shell.winfo_width(), 100)
                height = frame.winfo_reqheight() + 8
                shell.configure(height=height)
                backdrop.configure(height=height)
                backdrop.itemconfigure(window, width=width - 18)
                backdrop.delete("card")
                self.rounded_rect(backdrop, 1, 1, width - 1, height - 1, 17,
                                  fill=PANEL, outline="#28335f", width=1, tags="card")
                backdrop.tag_lower("card")

            frame.bind("<Configure>", fit)
            shell.bind("<Configure>", fit)
            self.label(frame, title.upper(), 10, CYAN, True).pack(anchor="w", pady=(0, 10))
            return frame

        def check(self, parent, name, variable, command, color=TEXT):
            box = tk.Checkbutton(parent, text=name, variable=variable, command=command,
                                 bg=PANEL, fg=color, selectcolor=CARD,
                                 activebackground=PANEL, activeforeground=TEXT,
                                 font=("Segoe UI", 10), highlightthickness=0,
                                 anchor="w", padx=0)
            return box

        def _build(self):
            header = tk.Frame(self.root, bg=BG, padx=24, pady=15)
            header.pack(fill="x")
            self.label(header, "POWER OVER FIBER", 21, TEXT, True).pack(side="left")
            self.label(header, "THERMAL LINK  /  INTERACTIVE SPECTRUM", 10, CYAN,
                       True).pack(side="left", padx=20, pady=7)
            identity = tk.Frame(header, bg=BG)
            identity.pack(side="right", padx=(16, 0))
            self.label(identity, "Nima Talebzadeh, Ph.D.", 12, "#b9a8ff", True).pack(anchor="e")
            self.label(identity, "POWER OVER FIBER", 8, MUTED, True).pack(anchor="e", pady=(2, 0))
            self.flux_header = self.label(header, "", 12, "#a990ff", True)
            self.flux_header.pack(side="right", padx=(0, 20))

            body = tk.Frame(self.root, bg=BG)
            body.pack(fill="both", expand=True, padx=18, pady=(0, 16))
            body.grid_columnconfigure(0, minsize=365)
            body.grid_columnconfigure(1, weight=1)
            body.grid_rowconfigure(0, weight=1)

            sidebar = tk.Frame(body, bg=BG, width=365)
            sidebar.grid(row=0, column=0, sticky="nsew", padx=(0, 14))
            sidebar.grid_propagate(False)
            scroll = tk.Canvas(sidebar, bg=BG, highlightthickness=0)
            bar = tk.Scrollbar(sidebar, command=scroll.yview)
            scroll.configure(yscrollcommand=bar.set)
            bar.pack(side="right", fill="y")
            scroll.pack(side="left", fill="both", expand=True)
            control = tk.Frame(scroll, bg=BG)
            window_id = scroll.create_window((0, 0), window=control, anchor="nw")
            control.bind("<Configure>", lambda e: scroll.configure(scrollregion=scroll.bbox("all")))
            scroll.bind("<Configure>", lambda e: scroll.itemconfigure(window_id, width=e.width))

            temp_frame = self.section(control, "Emitter temperature")
            self.temp_dial = tk.Canvas(temp_frame, height=74, bg=PANEL, highlightthickness=0)
            self.temp_dial.pack(fill="x")
            self.temp_dial.bind("<Configure>", lambda e: self.draw_temperature_dial())
            self.temp_dial.bind("<Button-1>", self.choose_temperature)
            self.draw_temperature_dial()

            distance_frame = self.section(control, "Fiber distance")
            self.distance_value = self.label(distance_frame, "1.00 km", 18, CYAN, True)
            self.distance_value.pack(anchor="w")
            tk.Scale(distance_frame, from_=0, to=15, resolution=0.01,
                     orient="horizontal", variable=self.distance, command=lambda _: self.schedule(),
                     showvalue=False, length=310, bg=PANEL, fg=TEXT,
                     troughcolor="#313d70", activebackground=CYAN,
                     highlightthickness=0, bd=0).pack(fill="x")
            ends = tk.Frame(distance_frame, bg=PANEL)
            ends.pack(fill="x")
            self.label(ends, "0 km", 9, MUTED).pack(side="left")
            self.label(ends, "15 km", 9, MUTED).pack(side="right")

            pv_frame = self.section(control, "PV band / wavelength range")
            self.check(pv_frame, "Select all PV cells", self.pv_all_var,
                       self.toggle_all_pv, CYAN).pack(anchor="w")
            for name, var in self.pv_vars.items():
                self.check(pv_frame, f"{name}    300–{PV_BANDS[name]} nm", var,
                           self.pv_selection_changed).pack(anchor="w")
            self.check(pv_frame, "Custom wavelength range", self.custom_var,
                       self.custom_changed, "#f596c8").pack(anchor="w", pady=(6, 0))
            custom_row = tk.Frame(pv_frame, bg=PANEL)
            custom_row.pack(fill="x", pady=(4, 0))
            self.low_entry = tk.Entry(custom_row, width=8, textvariable=self.low_var, bg=CARD,
                                      fg=TEXT, insertbackground=TEXT, relief="flat",
                                      font=("Segoe UI", 11), justify="center")
            self.low_entry.pack(side="left", ipady=5)
            self.label(custom_row, "to", 10, MUTED).pack(side="left", padx=8)
            self.high_entry = tk.Entry(custom_row, width=8, textvariable=self.high_var, bg=CARD,
                                       fg=TEXT, insertbackground=TEXT, relief="flat",
                                       font=("Segoe UI", 11), justify="center")
            self.high_entry.pack(side="left", ipady=5)
            self.label(custom_row, "nm", 10, MUTED).pack(side="left", padx=8)
            for entry in (self.low_entry, self.high_entry):
                entry.bind("<KeyRelease>", lambda e: self.schedule())
            self.custom_changed()

            fiber_frame = self.section(control, "Fiber type")
            self.fiber_all_var = tk.BooleanVar(value=True)
            self.check(fiber_frame, "Select all fibers", self.fiber_all_var,
                       self.toggle_all_fibers, CYAN).pack(anchor="w")
            for name, var in self.fiber_vars.items():
                self.check(fiber_frame, name, var, self.fiber_selection_changed,
                           FIBER_COLORS[name]).pack(anchor="w")

            photo_frame = self.section(control, "Power over fiber concept")
            self.photos = []
            for filename, caption in (("pof_scene.png", "Thermal source to remote receiver"),):
                path = ROOT / filename
                if path.exists():
                    img = Image.open(path).convert("RGB")
                    thumbnail = ImageOps.contain(img, (320, 145))
                    photo = ImageTk.PhotoImage(thumbnail)
                    self.photos.append(photo)
                    thumb = tk.Label(photo_frame, image=photo, bg=PANEL, cursor="hand2")
                    thumb.pack(pady=(4, 2))
                    thumb.bind("<Button-1>", lambda e, p=path: self.enlarge(p))
                    self.label(photo_frame, caption + "   •   click to enlarge", 9, MUTED).pack(pady=(0, 8))

            right = tk.Frame(body, bg=BG)
            right.grid(row=0, column=1, sticky="nsew")
            right.grid_rowconfigure(0, weight=1)
            right.grid_rowconfigure(1, weight=1)
            right.grid_columnconfigure(0, weight=1)

            plot_panel = tk.Frame(right, bg=PANEL, padx=10, pady=8)
            plot_panel.grid(row=0, column=0, sticky="nsew", pady=(0, 12))
            self.figure = Figure(figsize=(8, 4.6), dpi=100, facecolor=PANEL)
            self.ax = self.figure.add_subplot(111)
            self.chart = FigureCanvasTkAgg(self.figure, master=plot_panel)
            self.chart.get_tk_widget().pack(fill="both", expand=True)

            bottom = tk.Frame(right, bg=BG)
            bottom.grid(row=1, column=0, sticky="nsew")
            bottom.grid_columnconfigure(0, weight=3, uniform="bottom")
            bottom.grid_columnconfigure(1, weight=2, uniform="bottom")
            bottom.grid_rowconfigure(0, weight=1)

            result_panel = tk.Frame(bottom, bg=PANEL, padx=12, pady=9)
            result_panel.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
            self.label(result_panel, "TRANSMITTED IN-BAND OPTICAL POWER", 12, CYAN, True).pack(anchor="w")
            self.label(result_panel, "Ring: retention     •     Flux: W/cm²",
                       9, MUTED).pack(anchor="w", pady=(0, 4))
            result_canvas = tk.Canvas(result_panel, bg=PANEL, highlightthickness=0)
            result_scroll = tk.Scrollbar(result_panel, command=result_canvas.yview)
            result_canvas.configure(yscrollcommand=result_scroll.set)
            result_scroll.pack(side="right", fill="y")
            result_canvas.pack(fill="both", expand=True)
            self.result_content = tk.Frame(result_canvas, bg=PANEL)
            result_id = result_canvas.create_window((0, 0), window=self.result_content, anchor="nw")
            self.result_content.bind("<Configure>",
                                     lambda e: result_canvas.configure(scrollregion=result_canvas.bbox("all")))
            result_canvas.bind("<Configure>",
                               lambda e: result_canvas.itemconfigure(result_id, width=e.width))
            self.status = self.label(result_panel, "", 8, MUTED, wraplength=620, justify="left")
            self.status.pack(anchor="w", pady=(4, 0))

            schematic_panel = tk.Frame(bottom, bg=PANEL, padx=10, pady=9)
            schematic_panel.grid(row=0, column=1, sticky="nsew", padx=(8, 0))
            self.label(schematic_panel, "POWER OVER FIBER SYSTEM", 10, CYAN, True).pack(anchor="w")
            self.label(schematic_panel, "Emitter  →  optical cavity  →  fiber  →  PV", 8, MUTED).pack(
                anchor="w", pady=(3, 6))
            self.schematic_path = ROOT / "pof_schematic.png"
            self.schematic_source = Image.open(self.schematic_path).convert("RGB")
            self.schematic_view = tk.Label(schematic_panel, bg=PANEL, cursor="hand2")
            self.schematic_view.pack(fill="both", expand=True)
            self.schematic_view.bind("<Button-1>", lambda e: self.enlarge(self.schematic_path))
            self.schematic_panel = schematic_panel
            schematic_panel.bind("<Configure>", self.resize_schematic)

        def enlarge(self, path):
            from PIL import Image, ImageTk
            import tkinter as tk
            window = tk.Toplevel(self.root, bg=BG)
            window.title(path.stem)
            img = Image.open(path)
            img.thumbnail((1100, 750))
            photo = ImageTk.PhotoImage(img)
            view = tk.Label(window, image=photo, bg=BG)
            view.image = photo
            view.pack(padx=10, pady=10)

        def resize_schematic(self, event=None):
            if not hasattr(self, "schematic_source"):
                return
            from PIL import Image, ImageTk
            width = max(self.schematic_panel.winfo_width() - 24, 100)
            height = max(self.schematic_panel.winfo_height() - 58, 100)
            image = self.schematic_source.copy()
            image.thumbnail((width, height), Image.Resampling.LANCZOS)
            self.schematic_photo = ImageTk.PhotoImage(image)
            self.schematic_view.configure(image=self.schematic_photo)

        def toggle_all_pv(self):
            selected = self.pv_all_var.get()
            for var in self.pv_vars.values():
                var.set(selected)
            self.schedule()

        def pv_selection_changed(self):
            self.pv_all_var.set(all(var.get() for var in self.pv_vars.values()))
            self.schedule()

        def toggle_all_fibers(self):
            for var in self.fiber_vars.values():
                var.set(self.fiber_all_var.get())
            self.schedule()

        def fiber_selection_changed(self):
            self.fiber_all_var.set(all(var.get() for var in self.fiber_vars.values()))
            self.schedule()

        def custom_changed(self):
            state = "normal" if self.custom_var.get() else "disabled"
            for entry in (self.low_entry, self.high_entry):
                entry.configure(state=state)
            self.schedule()

        def schedule(self):
            if hasattr(self, "temp_dial"):
                self.draw_temperature_dial()
            if self.pending is not None:
                self.root.after_cancel(self.pending)
            self.pending = self.root.after(70, self.refresh)

        def draw_temperature_dial(self):
            if not hasattr(self, "temp_dial"):
                return
            canvas = self.temp_dial
            canvas.delete("all")
            width = max(canvas.winfo_width(), 300)
            values = list(TEMPERATURE_COLUMNS)
            centers = np.linspace(43, width - 43, len(values))
            for temp, center in zip(values, centers):
                selected = self.temperature.get() == temp
                color = CYAN if selected else "#4b5785"
                radius = 25
                canvas.create_oval(center-radius-4, 5, center+radius+4, 67,
                                   outline="#393b76" if selected else PANEL, width=2)
                canvas.create_oval(center-radius, 11, center+radius, 61,
                                   fill="#26305b" if selected else CARD,
                                   outline=color, width=2 if selected else 1)
                canvas.create_text(center, 32, text=f"{temp}", fill=TEXT,
                                   font=("Segoe UI", 10, "bold"))
                canvas.create_text(center, 48, text="K", fill=CYAN if selected else MUTED,
                                   font=("Segoe UI", 8, "bold"))

        def choose_temperature(self, event):
            width = max(self.temp_dial.winfo_width(), 300)
            centers = np.linspace(43, width - 43, len(TEMPERATURE_COLUMNS))
            index = min(range(len(centers)), key=lambda i: abs(centers[i] - event.x))
            self.temperature.set(list(TEMPERATURE_COLUMNS)[index])
            self.schedule()

        def selected_bands(self):
            bands = [(name, 300, cutoff) for name, cutoff in PV_BANDS.items()
                     if self.pv_vars[name].get()]
            if self.custom_var.get():
                try:
                    low, high = float(self.low_var.get()), float(self.high_var.get())
                except ValueError:
                    raise ValueError("Enter numeric custom wavelengths in nm")
                if not (300 <= low < high <= 2500) or high - low < 1:
                    raise ValueError("Custom range must be at least 1 nm within 300–2500 nm")
                bands.append(("Custom", low, high))
            return bands

        def draw_chart(self, temperature, distance, bands, fibers, input_curve, curves):
            ax = self.ax
            ax.clear()
            ax.set_facecolor(PANEL)
            for label, low, high in bands:
                if label == "Custom":
                    ax.axvspan(low, high, color=PV_COLORS[label], alpha=0.09)
                else:
                    ax.axvline(high, color=PV_COLORS[label], lw=1, ls=":", alpha=0.8)
                    ax.text(high + 8, 0.96, label, color=PV_COLORS[label], fontsize=8,
                            rotation=90, va="top")
            wl = model.wavelengths
            ax.plot(wl, input_curve, color=TEXT, linestyle="--", linewidth=1.6,
                    alpha=0.8, label="At fiber input")
            for name in fibers:
                ax.plot(wl, curves[name], color=FIBER_COLORS[name], linewidth=2.0,
                        label=name)
            ax.set(xlim=(300, 2500), ylim=(0, 1.05), xlabel="Wavelength (nm)",
                   ylabel="Spectral power / input peak",
                   title=f"Normalized spectral power  •  {temperature} K  •  {distance:.2f} km")
            ax.tick_params(colors=MUTED, labelsize=9)
            ax.xaxis.label.set_color(MUTED)
            ax.yaxis.label.set_color(MUTED)
            ax.title.set_color(TEXT)
            for spine in ax.spines.values():
                spine.set_color("#3d4877")
            ax.grid(color="#36416c", alpha=0.5)
            legend = ax.legend(facecolor=CARD, edgecolor="#465280", fontsize=8,
                               loc="upper right", ncol=min(2, len(fibers) + 1))
            for item in legend.get_texts():
                item.set_color(TEXT)
            self.figure.tight_layout(pad=1.6)
            self.chart.draw_idle()

        def draw_results(self, bands, fibers, results):
            import tkinter as tk
            for child in self.result_content.winfo_children():
                child.destroy()
            by_key = {(band, fiber): (fraction, flux) for band, fiber, fraction, flux in results}
            for name, low, high in bands:
                card_height = 48 + len(fibers) * 36
                card = tk.Canvas(self.result_content, height=card_height, bg=PANEL,
                                 highlightthickness=0, bd=0)
                card.pack(fill="x", pady=(0, 9))

                def paint(event, canvas=card, band_name=name, band_low=low,
                          band_high=high, selected_fibers=tuple(fibers), height=card_height):
                    width = max(event.width, 420)
                    canvas.delete("all")
                    self.rounded_rect(canvas, 5, 4, width - 5, height - 4, 18,
                                      fill=CARD, outline="#303c6b", width=1)
                    canvas.create_oval(16, 12, 32, 28, fill=PV_COLORS[band_name], outline="")
                    canvas.create_text(42, 20, text=band_name, fill=TEXT,
                                       font=("Segoe UI", 10, "bold"), anchor="w")
                    canvas.create_text(136, 20, text=f"{band_low:g}–{band_high:g} nm",
                                       fill=MUTED, font=("Segoe UI", 9), anchor="w")
                    canvas.create_text(width - 142, 20, text="RETAINED", fill=MUTED,
                                       font=("Segoe UI", 8, "bold"), anchor="e")
                    canvas.create_text(width - 22, 20, text="W/cm²", fill=MUTED,
                                       font=("Segoe UI", 8, "bold"), anchor="e")
                    for idx, fiber in enumerate(selected_fibers):
                        fraction, flux = by_key[band_name, fiber]
                        center_y = 48 + idx * 36
                        canvas.create_oval(16, center_y - 12, 40, center_y + 12,
                                           outline="#3a456d", width=4)
                        canvas.create_arc(16, center_y - 12, 40, center_y + 12,
                                          start=90, extent=-359.9 * min(max(fraction, 0), 1),
                                          style="arc", outline=FIBER_COLORS[fiber], width=4)
                        canvas.create_text(50, center_y, text=fiber, fill=TEXT,
                                           font=("Segoe UI", 9, "bold"), anchor="w")
                        canvas.create_text(width - 142, center_y,
                                           text=f"{100*fraction:.2f}%", fill=FIBER_COLORS[fiber],
                                           font=("Segoe UI", 13, "bold"), anchor="e")
                        canvas.create_text(width - 22, center_y, text=f"{flux:.4g}",
                                           fill=TEXT, font=("Segoe UI", 12, "bold"), anchor="e")
                        if idx < len(selected_fibers) - 1:
                            canvas.create_line(50, center_y + 17, width - 20, center_y + 17,
                                               fill="#354068", width=1)

                card.bind("<Configure>", paint)

        @staticmethod
        def rounded_rect(canvas, x1, y1, x2, y2, radius=16, **kwargs):
            r = min(radius, (x2 - x1) / 2, (y2 - y1) / 2)
            points = [x1+r, y1, x2-r, y1, x2, y1, x2, y1+r,
                      x2, y2-r, x2, y2, x2-r, y2, x1+r, y2,
                      x1, y2, x1, y2-r, x1, y1+r, x1, y1]
            return canvas.create_polygon(points, smooth=True, splinesteps=20, **kwargs)

        def refresh(self):
            self.pending = None
            temperature, distance = self.temperature.get(), self.distance.get()
            self.distance_value.configure(text=f"{distance:.2f} km")
            self.flux_header.configure(text=f"Emitter  {EMISSIVE_W_CM2[temperature]:.1f} W/cm²")
            try:
                bands = self.selected_bands()
                fibers = [name for name, var in self.fiber_vars.items() if var.get()]
                if not bands or not fibers:
                    raise ValueError("Select at least one PV band and one fiber")
                input_curve, curves, results = model.calculate(temperature, distance,
                                                               bands, fibers)
            except ValueError as exc:
                self.status.configure(text=str(exc), fg="#ff918d")
                return
            self.status.configure(
                text="W/cm² assumes all hemispherical emitter power is launched into the fiber. "
                     "Fiber coupling, PV absorptance, and electrical conversion are not included. "
                     "Blackbody tail >2400 nm and some fiber edges are extrapolated.",
                fg=MUTED)
            self.draw_chart(temperature, distance, bands, fibers, input_curve, curves)
            self.draw_results(bands, fibers, results)

        def run(self):
            self.root.mainloop()

    App().run()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Validate calculations without opening a window")
    args = parser.parse_args()
    model = SpectralModel()
    if args.check:
        bands = [(name, 300, cutoff) for name, cutoff in PV_BANDS.items()]
        for distance in (0.0, 1.0, 15.0):
            _, _, results = model.calculate(2000, distance, bands, list(FIBERS))
            print(f"Distance: {distance:g} km")
            for band, fiber, fraction, flux in results:
                print(f"  {band:9s} | {fiber:20s} | {100*fraction:8.4f}% | {flux:.6g} W/cm²")
        return
    start_gui(model)


if __name__ == "__main__":
    main()
