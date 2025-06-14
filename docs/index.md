# Canopy-App Documentation

<div class="hero-section">
  <img src="canopy-app-logo_no-bg.png" alt="Canopy-App Logo" style="max-height: 120px; margin-bottom: 1rem;">
  <p>Advanced atmospheric canopy modeling system for meteorology and air quality applications</p>
</div>

## Overview

The Canopy-App is a repository for low-level, stand-alone/column canopy parameterizations for testing and application to gridded atmospheric composition/air quality models. It provides detailed calculations for:

- **🌬️ Canopy meteorology** - Wind profiles, temperature, humidity, and wind adjustment factors
- **☀️ Radiation transfer** - Solar radiation attenuation and photolysis calculations
- **🌱 Biogenic emissions** - Volatile organic compound emissions from vegetation (MEGAN-based)
- **🍃 Dry deposition** - Gas and particle removal by vegetation surfaces
- **🌪️ Vertical diffusion** - In-canopy eddy diffusivities

**Authors:** Patrick Campbell, Zachary Moon, Wei-Ting Hung, Margaret Marvin, Quazi Rasool, and other NOAA research team members.

## Key Features

<div class="quick-links">
  <div class="quick-link-card">
    <h3>🔬 Comprehensive Physics</h3>
    <p>State-of-the-art parameterizations for canopy processes</p>
  </div>
  <div class="quick-link-card">
    <h3>📁 Multiple Formats</h3>
    <p>Supports both NetCDF and text input/output</p>
  </div>
  <div class="quick-link-card">
    <h3>⚙️ Flexible Configuration</h3>
    <p>Extensive namelist options for customization</p>
  </div>
  <div class="quick-link-card">
    <h3>📚 Well Documented</h3>
    <p>Complete API documentation and user guides</p>
  </div>
  <div class="quick-link-card">
    <h3>💻 Modern Fortran</h3>
    <p>Efficient Fortran 90+ implementation</p>
  </div>
  <div class="quick-link-card">
    <h3>🏆 Research Proven</h3>
    <p>Based on peer-reviewed scientific algorithms</p>
  </div>
</div>

## Quick Links

- [🚀 **Quick Start Guide**](getting-started/quickstart.md) - Get up and running quickly
- [📖 **Installation Guide**](getting-started/installation.md) - Detailed build instructions
- [⚙️ **Configuration Guide**](getting-started/configuration.md) - Setup and namelist options
- [🔬 **Science Documentation**](science/model-description.md) - Detailed model descriptions
- [💻 **API Reference**](api/overview.md) - Complete code documentation
- [📁 **Examples & Tutorials**](user-guide/examples.md) - Practical usage examples
- [📋 **Namelist Reference**](user-guide/namelist-reference.md) - Complete configuration options

## Current Components

The Canopy-App includes the following physics components:

### 1. **In-Canopy Winds and Wind Adjustment Factor (WAF)**
- **Purpose**: Wind profiles and wildfire spread applications
- **Namelist Options**: `ifcanwind`, `ifcanwaf`
- **Output Variables**: `ws` (m s⁻¹), `waf` (fraction)
- **Modules**: `canopy_wind_mod.F90`, `canopy_waf_mod.F90`
- **Reference**: Massman et al. (2017)

### 2. **In-Canopy Vertical Diffusion**
- **Purpose**: Eddy diffusivities for scaling resolved model layer diffusion
- **Namelist Option**: `ifcaneddy`
- **Output Variables**: `kz` (m² s⁻¹)
- **Module**: `canopy_eddyx_mod.F90`
- **Reference**: Massman et al. (2017), Makar et al. (2017)

### 3. **In-Canopy Photolysis Attenuation**
- **Purpose**: Scaling resolved model layer 1 photolysis rates
- **Namelist Option**: `ifcanphot`
- **Output Variables**: `rjcf` (fraction)
- **Module**: `canopy_phot_mod.F90`
- **Reference**: Massman et al. (2017), Makar et al. (2017)

### 4. **Leaf-Level Biogenic Emissions**
- **Purpose**: Volatile organic compound emissions (kg m⁻³ s⁻¹)
- **Namelist Option**: `ifcanbio`
- **Output Variables**: 19 species including isoprene, monoterpenes, sesquiterpenes
- **Module**: `canopy_bioemi_mod.F90`
- **Reference**: MEGANv2/v3 (Guenther et al., 2012), Clifton et al. (2021), Silva et al. (2020)

### 5. **Leaf-Level Gas Dry Deposition**
- **Purpose**: Gas removal by vegetation surfaces (cm s⁻¹)
- **Namelist Option**: `ifcanddepgas`
- **Output Variables**: RACM2 chemical mechanism species
- **Module**: `canopy_drydep_mod.F90`
- **Reference**: Zhang et al. (2003), ACCESS (Saylor 2013)

## Getting Started

### Prerequisites

**Required Dependencies:**
- Modern Fortran compiler (gfortran, ifort)
- NetCDF-Fortran Libraries (`-lnetcdf -lnetcdff`) for NetCDF I/O option

**Build Options:**
- `FC=gfortran` (default) or other compiler (e.g., `FC=ifort`)
- `DEBUG=0` (off, default), `DEBUG=1` (on), or `DEBUG=2` (extensive)
- `NC=0` (off) or `NC=1` (on, default) for NetCDF support

### Quick Build

```bash
# Basic build with debug flags and NetCDF
DEBUG=1 NC=1 make -C src

# Intel Fortran with debug flags
DEBUG=1 NC=1 FC=ifort make -C src

# Run the model
./canopy
```

### Configuration

Modify settings in the Fortran namelist file [`input/namelist.canopy`](https://github.com/noaa-oar-arl/canopy-app/blob/main/input/namelist.canopy) before running.

### Documentation Structure

1. **[Installation](getting-started/installation.md)** - Detailed build instructions and dependencies
2. **[Configuration](getting-started/configuration.md)** - Namelist options and input files
3. **[Quickstart](getting-started/quickstart.md)** - Get running in minutes
4. **[User Guide](user-guide/overview.md)** - Comprehensive usage documentation
5. **[Science](science/model-description.md)** - Physics and parameterizations
6. **[API Reference](api/overview.md)** - Complete code documentation
7. **[Examples](examples/basic.md)** - Practical examples and tutorials

## Recent Updates

!!! info "Latest Features"
    - Complete Doxygen documentation for all modules
    - Enhanced biogenic emission calculations with MEGAN v2/v3
    - Improved dry deposition parameterizations
    - Extended chemical mechanism support (RACM2)
    - Automated documentation deployment via GitHub Pages

## Documentation & Support

- **📚 [GitHub Repository](https://github.com/noaa-oar-arl/canopy-app)** - Source code and releases
- **🐛 [Report Issues](https://github.com/noaa-oar-arl/canopy-app/issues)** - Bug reports and feature requests
- **💬 [Discussions](https://github.com/noaa-oar-arl/canopy-app/discussions)** - Community discussions
- **🤝 [Contributing](development/contributing.md)** - Contribution guidelines

## Citation

If you use Canopy-App in your research, please cite:

```bibtex
@software{canopy_app_2024,
  title={Canopy-App: Atmospheric Canopy Modeling System},
  author={Campbell, P.C. and Moon, Z. and Hung, W.T. and Marvin, M. and Rasool, Q. and others},
  year={2024},
  url={https://github.com/noaa-oar-arl/canopy-app},
  doi={10.5281/zenodo.8403649}
}
```

---

*Developed by NOAA Air Resources Laboratory and contributors*
