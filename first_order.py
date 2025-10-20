import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
muon_count = 30000      # number of simulated muons
energy_min = 0.5        # GeV, lower energy bound for sampling
energy_max = 20.0       # GeV, upper energy bound
k = 0.0136              # 13.6 MeV converted to GeV for scattering formula (Highland formula constant)

# Acceptance threshold: angular cutoff for muon acceptance
# 0.03 rad ≈ 1.7° is typical for cosmic-ray muon detectors
# Rationale:
#   - Small enough to reject most muons scattered by dense materials (iron)
#   - Large enough to accept most muons through low-density materials (ice)
#   - Realistic for scintillator/tracker detector angular resolution (~0.1-0.5 mrad precision)
#   - Creates clear tomographic contrast: ~50% reduction between ice and iron
# Reference: Standard choice in muon tomography experiments (e.g., MINOS, KamLAND)
acceptance_threshold = 0.03  # radians (~1.7 degrees)

# Material properties: radiation lengths (g/cm^2) and densities (g/cm^3)
# Source: PDG 2024 (https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf)
X0_ice = 36.1     # Radiation length of ice (H₂O)
X0_iron = 13.84   # Radiation length of iron (Fe)
rho_ice = 0.92    # Density of ice at 0°C (g/cm³)
rho_iron = 7.87   # Density of iron at room temperature (g/cm³)

# Geometry of the phantom: 20 cm of ice, with 16 cm of iron at the centre
thickness_ice_total = 20.0
thickness_block = 16.0
thickness_ice_above = (thickness_ice_total - thickness_block) / 2.0
thickness_ice_below = thickness_ice_above

# Convert material thicknesses to g/cm^2 and compute x/X0 for each configuration
x_ice = thickness_ice_total * rho_ice
x_over_X0_ice = x_ice / X0_ice
# Combined material thickness for iron+ice (weighted by radiation length)
x_over_X0_combined = (
    thickness_ice_above * rho_ice / X0_ice +
    thickness_block * rho_iron / X0_iron +
    thickness_ice_below * rho_ice / X0_ice
)

# Helper function to sample muon energies from an E^-2.7 distribution
# Cosmic muon spectrum at sea level: dN/dE ∝ E^-2.7 (Gaisser 2004, PDG 2024)
# See also: https://indico.cern.ch/event/975141/contributions/4137563/attachments/2156087/3646193/2020-Lecture-3-Interactions%20of%20Particles%20with%20Matter.pdf (2.3.2 Cosmic Ray Muon Spectrum)
emin_pow = energy_min**(-1.7)
emax_pow = energy_max**(-1.7)
def sample_muon_energies(n):
    """Sample n muon energies from E^-2.7 cosmic spectrum using inverse transform method."""
    u = np.random.rand(n)
    return (emin_pow + (emax_pow - emin_pow) * u)**(-1.0 / 1.7)

# Sampling for scattering angle distributions
# This is a first-order approximation, not full Monte Carlo.
#   - Scattering angles are sampled randomly (2D Gaussian) but energy loss is not sampled (deterministic via Highland formula)
#   - Material traversal is not simulated (1D path assumed, no angular deflection)
#   - Muon trajectories do not evolve (straight line through material)
# For full Monte Carlo: use GEANT4 (tracks muon step-by-step through detector)
def compute_acceptance(threshold, use_block=False, count=muon_count):
    energies = sample_muon_energies(count)
    momenta = energies  # assuming beta ≈ 1
    if use_block:
        x_ratio = x_over_X0_combined
    else:
        x_ratio = x_over_X0_ice
    theta0 = k / momenta * np.sqrt(x_ratio)
    # Draw scattering angles in x and y (2D Gaussian). compute radial scattering angle
    gauss = np.random.randn(count, 2)
    scattering = np.sqrt((theta0[:, None] * gauss)**2).sum(axis=1)
    return np.mean(scattering < threshold)

# Compute acceptance for pure ice and iron+ice
acc_ice = compute_acceptance(acceptance_threshold, use_block=False)
acc_block = compute_acceptance(acceptance_threshold, use_block=True)

# Define a scanning grid: 20×20 cm phantom divided into 5‑cm cells
# Grid cells are centered at: [-7.5, -2.5, 2.5, 7.5] (4×4 grid, 5 cm spacing)
grid_positions = [-7.5, -2.5, 2.5, 7.5]
acceptance_map = np.ones((4, 4))
relative_reduction_map = np.zeros_like(acceptance_map)

# Populate the map: only the center cell (1 out of 4×4) contains the iron block
# Iron block: 5×5 cm centered at origin, bounded by ±2.5 cm (strict inequality)
for ix, x in enumerate(grid_positions):
    for iy, y in enumerate(grid_positions):
        if abs(x) < 2.5 and abs(y) < 2.5:
            # Muons traverse the iron block here (center cell only)
            acceptance_map[iy, ix] = acc_block
        else:
            acceptance_map[iy, ix] = acc_ice
relative_reduction_map = 1.0 - (acceptance_map / acc_ice)

# Plot the acceptance reduction heatmap. Styles could still be improved.
plt.figure(figsize=(8, 6.5))
plt.style.use('seaborn-v0_8-darkgrid')
extent = [grid_positions[0] - 2.5, grid_positions[-1] + 2.5,
          grid_positions[0] - 2.5, grid_positions[-1] + 2.5]
im = plt.imshow(relative_reduction_map, extent=extent, origin='lower',
                cmap='YlOrRd', vmin=0, vmax=np.max(relative_reduction_map),
                interpolation='nearest', alpha=0.95)
cbar = plt.colorbar(im, label='Relative reduction in acceptance', pad=0.02)
cbar.ax.tick_params(labelsize=10)
plt.title('Simulated reduction in muon acceptance due to iron block', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Scan position x (cm)', fontsize=12, fontweight='bold')
plt.ylabel('Scan position y (cm)', fontsize=12, fontweight='bold')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('acceptance_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()

# Generate and plot scattering-angle distributions for illustration.
energies = sample_muon_energies(50000)
momenta = energies
theta0_ice = k / momenta * np.sqrt(x_over_X0_ice)
theta0_block = k / momenta * np.sqrt(x_over_X0_combined)
scatt_ice = np.sqrt((theta0_ice[:, None] * np.random.randn(50000, 2))**2).sum(axis=1)
scatt_block = np.sqrt((theta0_block[:, None] * np.random.randn(50000, 2))**2).sum(axis=1)

plt.figure(figsize=(10, 6))
plt.style.use('seaborn-v0_8-darkgrid')
bins = np.linspace(0, 0.1, 50)
plt.hist(scatt_ice, bins=bins, alpha=0.6, label='Pure ice', color='#3498db', edgecolor='black', linewidth=1.2)
plt.hist(scatt_block, bins=bins, alpha=0.6, label='Iron inclusion', color='#e74c3c', edgecolor='black', linewidth=1.2)
plt.axvline(acceptance_threshold, color='#2ecc71', linestyle='--', linewidth=2.5,
            label=f'Acceptance threshold ({acceptance_threshold:.3f} rad)')
plt.xlabel('Scattering angle (rad)', fontsize=12, fontweight='bold')
plt.ylabel('Counts', fontsize=12, fontweight='bold')
plt.title('Distribution of simulated scattering angles', fontsize=14, fontweight='bold', pad=15)
plt.legend(fontsize=11, loc='upper right', framealpha=0.95, edgecolor='black', fancybox=True, shadow=True)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
plt.tight_layout()
plt.savefig('scattering_distribution.png', dpi=150, bbox_inches='tight')
plt.close()

print(f'Pure-ice acceptance: {acc_ice:.3f}') # Fraction of muons with scattering angle < 0.03 rad after traversing 20 cm ice
print(f'Iron+ice acceptance: {acc_block:.3f}') # Fraction of muons with scattering angle < 0.03 rad after traversing iron+ice
print(f'Relative reduction: {1 - acc_block / acc_ice:.3f}') # Relative reduction in acceptance due to iron block   


