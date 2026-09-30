<img align="right" width="230" src="docs/images/pof_hero.png" alt="Thermal source transferring power to drones over fiber">

# Thermal Power over Fiber Dashboard

Interactive desktop tool for estimating how much blackbody emission remains in selected photovoltaic wavelength bands after transmission through an optical fiber.

**Nima Talebzadeh, Ph.D.**  
Python · Tkinter · NumPy · Matplotlib

<br clear="right">

## Dashboard preview

![Power over fiber dashboard preview](docs/images/dashboard_preview.png)

## What it does

- Select a blackbody temperature: 1500, 2000, or 2500 K.
- Set fiber distance from 0 to 15 km.
- Select Si, GaSb, InGaAsSb, any combination, or a custom wavelength interval.
- Select one fiber or compare the Purchased, UHV, and low attenuation (OH) fibers.
- View the input and transmitted spectra on one normalized scale.
- Compare in-band optical retention and estimated transmitted spectral exitance in W/cm².

## Physical model

Attenuation \(\alpha(\lambda)\) is in dB/km. At wavelength \(\lambda\) and distance \(L\) in km, the model calculates transmission as

\[
T(\lambda,L)=10^{-\alpha(\lambda)L/10}.
\]

The in-band retention shown in each result is

\[
R_{\mathrm{band}}=\frac{\int_{\lambda_{\min}}^{\lambda_{\max}}S(\lambda)T(\lambda,L)\,d\lambda}
{\int_{\lambda_{\min}}^{\lambda_{\max}}S(\lambda)\,d\lambda},
\]

where \(S(\lambda)\) is the input blackbody spectrum. The default bands are 300–1100 nm for Si, 300–1700 nm for GaSb, and 300–2500 nm for InGaAsSb. A custom interval can be selected within 300–2500 nm.

The displayed W/cm² estimate weights the transmitted spectrum by Planck's law and scales it to the supplied total emissive exitances: 28.7 W/cm² at 1500 K, 90.7 W/cm² at 2000 K, and 221.5 W/cm² at 2500 K. It assumes all hemispherical emission is launched into the fiber. It does not include fiber coupling loss, receiver absorptance, PV temperature effects, or electrical conversion efficiency. It represents an optical flux per emitter area, not the actual power in watts at the fiber output.

## Data handling and limits

- Blackbody spectra are read from `blackbody.xlsx`; the supplied spectra match Planck's law at their measured wavelengths. The 2400–2500 nm tail follows Planck's law.
- Fiber attenuation data are read from `attenuation.xlsx`, interpolated between measurements, and linearly extrapolated beyond measured wavelength limits. Negative extrapolated attenuation is clipped to zero.
- Attenuation is smoothed with a Gaussian kernel with 20 nm standard deviation before plotting and integration. This reduces narrow measured variations and can soften real absorption features.
- Results that depend on wavelengths outside measured fiber data are estimates. Dashed segments in the attenuation plot mark extrapolated regions.
- The spectrum plot is normalized to the input spectrum's peak so transmitted curves can be compared directly with the input curve.

## Download the Windows application

The repository builds a Windows `.exe` automatically with GitHub Actions. To download it:

1. Open the repository's **Actions** tab and select the latest **Build Windows app** run.
2. Download the `PowerOverFiberDashboard-Windows` artifact.
3. Extract the downloaded artifact ZIP and run `PowerOverFiberDashboard.exe`.

When a version tag such as `v1.0.0` is pushed, the workflow also attaches the executable to a GitHub Release.

To publish the first release after pushing the source repository:

```powershell
git tag v1.0.0
git push origin v1.0.0
```

GitHub Actions builds the Windows `.exe` and adds it to the release for that tag. Others can download it from the repository's **Releases** page.

## Build or run from source

Requires Python 3.12 or newer on Windows. From this folder:

```powershell
python -m pip install -r requirements.txt
python fiber_thermal_gui.py
```

To build the one-file Windows executable locally:

```powershell
python -m pip install -r requirements-build.txt
./build_windows.ps1
```

The executable is created at `dist/PowerOverFiberDashboard.exe`. The workbook data and image assets are bundled into it; users do not need a separate Python installation.

## Check the calculations

Run the model checks at 0, 1, and 15 km without opening the GUI:

```powershell
python fiber_thermal_gui.py --check
```

At zero distance, each in-band retention should be 100%. The values then decrease with distance according to the selected fiber's attenuation spectrum.

## Repository contents

- `fiber_thermal_gui.py` — desktop interface and spectral model
- `blackbody.xlsx`, `attenuation.xlsx` — source spectra and fiber attenuation data
- `pof_scene.png`, `pof_schematic.png` — interface illustrations
- `.github/workflows/build-windows.yml` — Windows build and release workflow
