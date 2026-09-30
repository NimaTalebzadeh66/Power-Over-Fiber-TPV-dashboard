$ErrorActionPreference = "Stop"

python -m PyInstaller `
    --noconfirm `
    --clean `
    --onefile `
    --windowed `
    --name PowerOverFiberDashboard `
    --icon pof_app.ico `
    --add-data "blackbody.xlsx;." `
    --add-data "attenuation.xlsx;." `
    --add-data "pof_scene.png;." `
    --add-data "pof_schematic.png;." `
    --collect-data matplotlib `
    --hidden-import matplotlib.backends.backend_tkagg `
    fiber_thermal_gui.py

Write-Host "Built dist\PowerOverFiberDashboard.exe"
