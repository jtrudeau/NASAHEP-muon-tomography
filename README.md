# Muon Tomography Simulation

A first-order Monte Carlo simulation that demonstrates muon tomography – using cosmic-ray muons to image subsurface density by measuring scattering.

## How it works roughly:

1. Cosmic muons traverse ice and or an iron block
2. Iron scatters muons more than ice (4× stronger deflection)
3. Detector counts muons with small scattering angles (acceptance threshold estimated at 0.03 rad)
4. Central region (iron) has lower acceptance which creates contrast map
5. Map reveals the location & shape of the dense iron block

## The Physics

- Muons undergo multiple Coulomb scattering in matter
- Highland formula predicts scattering angles: `θ₀ = (k/p) √(x/X₀)`
- Radiation length X₀ determines scattering power:
  - Ice: X₀ = 36.1 g/cm² (low density)
  - Iron: X₀ = 13.84 g/cm² (higher density)


## Outputs of the simulation

- `acceptance_heatmap.png`– 2D map showing relative acceptance reduction (bright = iron, dark = ice)
- `scattering_distribution.png` – Histograms comparing scattering angles (pure ice vs. iron+ice)

## Expected Results

- Pure ice acceptance: ~85%
- Iron+ice acceptance: ~40%
- Relative reduction: ~50% (for clear tomographic contrast)

## Real-world applications?

- Ice moon subsurface imaging (Europa, Enceladus)
- Geological surveys and archaeological scanning
- Mars...


Takeaway at first order: Denser materials scatter muons more. By measuring the scattering signature across a spatial grid, we can map subsurface density without drilling.
