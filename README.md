# viscoprony
**viscoprony** is a Python package for viscoelastic material modeling. It generates Prony-series coefficients from dynamic mechanical analysis (DMA) data to describe time-dependent behavior of polymers and other viscoelastic materials.

# Purpose
This tool streamlines the process of polymer and material characterization for FEA implementations (e.g., ABAQUS, ANSYS)

1. Building master curves using Time-Temperature Superposition (TTS)
2. Fitting Prony-series parameters (relaxation moduli and times) from DMA data
3. Visualizing shift factors, storage/loss moduli, and fitted curves

# Clone the repository
git clone [https://github.com/yourusername/viscoprony.git](https://github.com/MrJSimon/viscoprony.git)

# Navigate into the project folder
cd viscoprony

# Install in editable/development mode
pip install -e .

# Visualizations

**Shift Factors vs Temperature**: Logarithmic shift factors $a_T$ calculated from DMA data. 
<p align="center">
  <img src="./docs/images/ShiftfactorVsTemperatures.png" alt="Shift Factors" width="70%">
</p>

**Storage and Loss Moduli**: Raw DMA curves showing $E'$ (storage modulus) and $E''$  (loss modulus).
<p align="center">
  <img src="./docs/images/StorageLossCurves.png" alt="Storage and loss moduli curves from DMA data" width="70%">
</p>

**Master Curve**: Constructed master curve showing the material’s behavior over a wide frequency range.
<p align="center">
  <img src="./docs/images/MasterCurve.png" alt="Master curve of viscoelastic material" width" width="70%">
</p>

**Prony-Series Fit**: Fitted Prony-series overlaying the master curve for validation.
<p align="center">
  <img src="./docs/images/PronyseriesFit.png" alt="Prony-series fit on master curve" width" width="70%">
</p>
