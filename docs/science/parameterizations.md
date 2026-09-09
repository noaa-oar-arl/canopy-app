# Parameterizations

The Canopy-App model includes various parameterizations for physical and chemical processes within and above forest canopies. This section provides detailed descriptions of the scientific formulations and parameters used.

## Turbulence Parameterization

### K-Theory Approach

The model uses K-theory to parameterize turbulent transport:

$$
\overline{w'\phi'} = -K_{\phi} \frac{\partial \phi}{\partial z}
$$

where:
- $\overline{w'\phi'}$ is the turbulent flux of scalar $\phi$
- $K_{\phi}$ is the eddy diffusivity for scalar $\phi$
- $z$ is height above ground

### Mixing Length Model

The eddy diffusivity is calculated using a mixing length approach:

$$
K_m = l_m^2 \sqrt{\left(\frac{\partial U}{\partial z}\right)^2 + \left(\frac{\partial V}{\partial z}\right)^2}
$$

$$
K_h = \frac{K_m}{\Pr_t}
$$

where:
- $K_m$ is momentum diffusivity
- $K_h$ is heat/scalar diffusivity
- $l_m$ is mixing length
- $\Pr_t$ is turbulent Prandtl number (typically 0.7-1.3)

### Mixing Length Formulation

Within the canopy, the mixing length is parameterized as:

$$
l_m(z) = \begin{cases}
\beta h \left(\frac{z}{h}\right)^n & \text{for } z < h \\
\kappa z & \text{for } z \geq h
\end{cases}
$$

where:
- $h$ is canopy height
- $\beta$ is scaling parameter (typically 0.1-0.3)
- $n$ is shape parameter (typically 1-3)
- $\kappa$ is von Karman constant (0.41)

## Radiation Parameterization

<!-- ### Two-Stream Approximation

Solar radiation transfer uses the two-stream approximation:

$$
\mu \frac{dI^+}{d\tau} = I^+ - \omega \beta_0 I^+ - \omega \beta_1 I^-
$$

$$
-\mu \frac{dI^-}{d\tau} = I^- - \omega \beta_1 I^+ - \omega \beta_0 I^-
$$

where:
- $I^+, I^-$ are upward and downward radiation intensities
- $\tau$ is optical depth
- $\omega$ is single scattering albedo
- $\beta_0, \beta_1$ are phase function parameters
- $\mu = \cos(\theta)$ for solar zenith angle $\theta$ -->

### Leaf Area Distribution

The cumulative leaf area index from canopy top:

$$
LAI(z) = LAI_{total} \exp\left(-\alpha \left(\frac{h-z}{h}\right)^{\beta}\right)
$$

Parameters:
- $\alpha$: extinction coefficient (typically 0.5-2.0)
- $\beta$: shape parameter (typically 0.5-2.0)

## Biogenic Emissions

### Simplified Parameterization

This compact workflow summarizes the scientific parameterization and principal output
choices for use as a one-page or journal-article figure.

```mermaid
flowchart LR
	classDef input fill:#FFF3CD,stroke:#B77900,color:#3D2900,stroke-width:2px
	classDef model fill:#E3F2E8,stroke:#287A45,color:#123A22,stroke-width:2px
	classDef factor fill:#F8E5EF,stroke:#A93669,color:#4D1831,stroke-width:2px
	classDef decision fill:#EEE7F8,stroke:#68429A,color:#2F1948,stroke-width:2px
	classDef output fill:#DFF3F5,stroke:#087987,color:#073B41,stroke-width:2px

	subgraph A["Inputs and initialization"]
		direction TB
		CFG["Biogenic options<br/>species, vertical mode,<br/>gamma and loss switches"]:::input
		STATE["Canopy and environment<br/>LAI profile, vegetation, PPFD,<br/>Tleaf, soil, CO2, O3, wind"]:::input
		HIST["Instantaneous or<br/>24 h / 240 h history"]:::input
	end

	subgraph B["Species activity"]
		direction TB
		PARM["Species x vegetation parameters<br/>emission factor, light fraction,<br/>temperature and stress coefficients"]:::model
		ENV["Canopy environment<br/>gamma(T, PPFD)<br/>sunlit + shaded leaves"]:::factor
		MOD["Response factors<br/>gamma(CO2) x gamma(leaf age)<br/>x gamma(soil moisture)<br/>x gamma(AQ, heat, cold, wind)"]:::factor
		EMIS["Potential emission<br/>EF x gamma(T, PPFD) x response factors<br/>x canopy environment coefficient"]:::model
	end

	subgraph C["Vertical treatment"]
		direction TB
		VERT{"biovert_opt"}:::decision
		PROFILE["0 | Layer-resolved source<br/>weighted by leaf area density"]:::model
		INTEGRATE["1-3 | Canopy-integrated flux<br/>plant, Gaussian, or equal weighting<br/>x eligible canopy loss factor"]:::model
	end

	subgraph D["Model output"]
		direction TB
		VOL["3D volumetric emissions<br/>kg m-3 s-1"]:::output
		AREA["2D areal flux in top layer<br/>kg m-2 s-1"]:::output
		WRITE["1D point text or<br/>2D gridded NetCDF"]:::output
	end

	CFG --> PARM
	STATE --> ENV
	HIST --> ENV
	PARM --> ENV
	PARM --> MOD
	ENV --> EMIS
	MOD --> EMIS
	EMIS --> VERT
	VERT -->|0| PROFILE --> VOL --> WRITE
	VERT -->|1-3| INTEGRATE --> AREA --> WRITE

	linkStyle default stroke:#4F6470,stroke-width:1.8px
```

**Figure 1.** Simplified Canopy-App biogenic emissions parameterization. Environmental
activity and species-specific response factors modify the mapped emission factor before
the source is retained as a layer-resolved profile or vertically integrated to a
top-canopy flux.

### Parameterization Workflow

The workflow below follows the implemented path from namelist initialization through
[`canopy_calcs.F90`](../../src/canopy_calcs.F90),
[`canopy_bioemi_mod.F90`](../../src/canopy_bioemi_mod.F90), and the text or NetCDF
writers. Colors group configuration, environmental inputs, activity factors, vertical
treatment, unit conversion, and output.

```mermaid
flowchart TB
	classDef config fill:#E8F1FA,stroke:#1B5E8C,color:#123047,stroke-width:2px
	classDef input fill:#FFF4D6,stroke:#C47A00,color:#4A3000,stroke-width:2px
	classDef process fill:#E7F6EC,stroke:#287A45,color:#123A22,stroke-width:2px
	classDef gamma fill:#FCE8F1,stroke:#B23A6F,color:#541A35,stroke-width:2px
	classDef decision fill:#F3E8FF,stroke:#7048A8,color:#321C50,stroke-width:2px
	classDef units fill:#FFF0E5,stroke:#C45118,color:#54230C,stroke-width:2px
	classDef output fill:#E5F6F8,stroke:#087D8D,color:#073C43,stroke-width:2px
	classDef stop fill:#FDE8E7,stroke:#B8322A,color:#551814,stroke-width:2px

	subgraph INIT["1 | Configuration and initialization"]
		direction LR
		NML["canopy_readnml<br/>ifcanbio, biospec_opt<br/>biovert_opt, bio_cce<br/>gamma switches, loss_opt"]:::config
		ALLOC["canopy_alloc<br/>allocate selected species arrays"]:::config
		INITARR["canopy_init<br/>initialize instantaneous, 24 h,<br/>240 h, stress, and emission arrays"]:::config
		NML --> ALLOC --> INITARR
	end

	subgraph CALCS["2 | Driver preparation in canopy_calcs.F90"]
		direction LR
		MET["Canopy state and meteorology<br/>LAI, FCLAI, canopy height, Tair,<br/>soil moisture, CO2, u*, O3, wind"]:::input
		RAD["Sunlit / shaded leaf environment<br/>FSUN, PPFD, Tleaf"]:::input
		HIST{"hist_opt"}:::decision
		INST["0: use instantaneous<br/>PPFD, Tleaf, T2, wind"]:::process
		AVG["1: update 24 h / 240 h means<br/>and daily stress extrema"]:::process
		LEAF["leafage_opt=0<br/>update past/current LAI<br/>at lai_tstep"]:::process
		SELECT{"ifcanbio and<br/>vegetated canopy?"}:::decision
		SPEC["biospec_opt<br/>0 = all 19 species<br/>1...19 = selected species"]:::config

		MET --> HIST
		RAD --> HIST
		HIST -->|0| INST
		HIST -->|1| AVG
		INST --> SELECT
		AVG --> SELECT
		LEAF --> SELECT
		SPEC --> SELECT
	end

	INITARR --> MET
	SELECT -->|No| ZERO["Set selected emission<br/>profiles to zero"]:::stop
	SELECT -->|Yes| BIO["CANOPY_BIO for each selected species<br/>EMI_IND = 1...19"]:::process

	subgraph PARAM["3 | Species parameterization in canopy_bioemi_mod.F90"]
		direction TB
		BIOP["CANOPY_BIOP<br/>Map EMI_IND + land use + vegetation to<br/>EF, LDF, BETA, CT1, CEO and response coefficients"]:::process
		TL["Temperature activity<br/>Eopt and Topt from 24 h / 240 h Tleaf<br/>sunlit + shaded LDF and LIF terms"]:::gamma
		LIGHT["Light activity<br/>alpha and Cp from 24 h / 240 h PPFD<br/>sunlit + shaded PPFD response"]:::gamma
		ENV["Canopy environment activity<br/>gamma_T,PPFD = LDF gamma_LDF<br/>+ (1-LDF) gamma_LIF"]:::gamma
		GCO2["gamma_CO2<br/>isoprene only; otherwise 1"]:::gamma
		GLEAF["gamma_leafage<br/>new + growing + mature + old foliage"]:::gamma
		GSOIL["gamma_soil moisture<br/>root-weighted soil layers"]:::gamma
		GSTRESS["Stress gammas<br/>gamma_AQ x gamma_HT x gamma_LT x gamma_HW"]:::gamma
		PRODUCT["Layer activity product<br/>EF x gamma_T,PPFD x gamma_CO2 x CCE<br/>x gamma_leafage x gamma_soil<br/>x gamma_AQ x gamma_HT x gamma_LT x gamma_HW"]:::process

		BIO --> BIOP
		BIOP --> TL
		BIOP --> LIGHT
		TL --> ENV
		LIGHT --> ENV
		BIOP --> GCO2
		BIOP --> GLEAF
		BIOP --> GSOIL
		BIOP --> GSTRESS
		ENV --> PRODUCT
		GCO2 --> PRODUCT
		GLEAF --> PRODUCT
		GSOIL --> PRODUCT
		GSTRESS --> PRODUCT
	end

	subgraph VERTICAL["4 | Vertical treatment, canopy loss, and conversion"]
		direction TB
		VERT{"biovert_opt"}:::decision
		V0["0 | Full 3D profile<br/>FLAI = delta FCLAI x LAI / dz<br/>E(z) = FLAI(z) x activity product"]:::process
		V1["1 | Plant-distribution weighting<br/>VPGWT = FLAI / sum(FLAI)"]:::process
		V2["2 | MEGAN-like five-level<br/>Gaussian weighting"]:::process
		V3["3 | Equal weighting<br/>VPGWT = 1 / canopy layers"]:::process
		LOSSSEL{"Integrated path and<br/>LOSSIND = 0 or EMI_IND?"}:::decision
		LOSS0["loss_opt=0<br/>CANLOSS = 1"]:::process
		LOSS1["loss_opt=1<br/>lifetime / u* / canopy-height loss"]:::process
		LOSS2["loss_opt=2<br/>CANLOSS = loss_set"]:::process
		INT["Top-canopy areal flux<br/>LAI x EF x weighted gamma product<br/>x CANLOSS"]:::process
		U3["Multiply by 2.7777778e-13<br/>microgram m-3 h-1 to kg m-3 s-1"]:::units
		U2["Multiply by 2.7777778e-13<br/>microgram m-2 h-1 to kg m-2 s-1<br/>store in top model layer"]:::units

		PRODUCT --> VERT
		VERT -->|0| V0 --> U3
		VERT -->|1| V1 --> LOSSSEL
		VERT -->|2| V2 --> LOSSSEL
		VERT -->|3| V3 --> LOSSSEL
		LOSSSEL -->|No| LOSS0
		LOSSSEL -->|Yes, option 0| LOSS0
		LOSSSEL -->|Yes, option 1| LOSS1
		LOSSSEL -->|Yes, option 2| LOSS2
		LOSS0 --> INT
		LOSS1 --> INT
		LOSS2 --> INT
		INT --> U2
	end

	subgraph OUTPUT["5 | Output pathway"]
		direction LR
		OUTTYPE{"Build / input mode"}:::decision
		TXT["1D point-list pathway<br/>canopy_write_txt<br/>*_bio.txt; all species required<br/>lat, lon, z, LAD, emission profiles"]:::output
		NCF["2D gridded pathway<br/>canopy_write_ncf<br/>NetCDF fields on lon x lat x layer x time"]:::output
		OUTTYPE -->|Text build| TXT
		OUTTYPE -->|NETCDF build| NCF
	end

	U3 --> OUTTYPE
	U2 --> OUTTYPE

	linkStyle default stroke:#536878,stroke-width:1.6px
```

**Figure 2.** Detailed Canopy-App biogenic emissions implementation workflow. The parameterization is
evaluated independently for each selected species. `biospec_opt=0` evaluates all 19
species; values 1--19 select an individual species.

!!! note "Vertical option and units"
    `biovert_opt=0` retains the layer-resolved source and converts
    $\text{microgram m}^{-3}\text{ h}^{-1}$ to $\text{kg m}^{-3}\text{ s}^{-1}$.
    Options 1--3 vertically integrate the source, apply the eligible canopy loss
    factor, convert $\text{microgram m}^{-2}\text{ h}^{-1}$ to
    $\text{kg m}^{-2}\text{ s}^{-1}$, and place the result in the top model layer.
    Current text and NetCDF variable metadata label biogenic fields as volumetric
    emissions, so users should interpret integrated-option output according to
    `biovert_opt`.

Editable Mermaid sources are available for the
[simplified journal figure](../development/biogenic_emissions_flowchart_simplified.mmd)
and the [detailed implementation flowchart](../development/biogenic_emissions_flowchart.mmd)
for vector export to SVG or PDF.

### Emission Calculation

Following Guenther et al. (2012):

$$
E_i = \epsilon_i \cdot \gamma_{T,P} \cdot \gamma_{CO_2} \cdot C_{CE}
\cdot \gamma_{age} \cdot \gamma_{SM} \cdot \gamma_{AQ}
\cdot \gamma_{HT} \cdot \gamma_{LT} \cdot \gamma_{HW}
$$

where:
- $\epsilon_i$: species- and vegetation-dependent emission factor
- $\gamma_{T,P}$: combined temperature and PPFD activity factor
- $\gamma_{CO_2}$: carbon dioxide inhibition factor for isoprene (unity otherwise)
- $C_{CE}$: canopy environment coefficient (`bio_cce`)
- $\gamma_{age}$: leaf age activity factor
- $\gamma_{SM}$: soil moisture activity factor
- $\gamma_{AQ}$: air-quality stress factor
- $\gamma_{HT}$ and $\gamma_{LT}$: high- and low-temperature stress factors
- $\gamma_{HW}$: high-wind stress factor

### Temperature and Light Dependence

$$
\gamma_{T,P} = C_T \cdot C_L
$$

$$
C_T = \frac{\exp\left(\frac{C_{T1}(T-T_s)}{RT_sT}\right)}{1 + \exp\left(\frac{C_{T2}(T-T_{M})}{RT_sT}\right)}
$$

$$
C_L = \frac{\alpha C_L1 PPFD}{\sqrt{1 + \alpha^2 PPFD^2}}
$$

Parameters:
- $C_{T1} = 95,000$ J/mol
- $C_{T2} = 230,000$ J/mol
- $T_s = 303$ K (standard temperature)
- $T_M = 314$ K (maximum temperature)
- $\alpha = 0.0027$ (empirical coefficient)
- $C_{L1} = 1.066$ (empirical coefficient)

## Dry Deposition

### Resistance Model

Total deposition velocity:

$$
v_d = \frac{1}{r_a + r_b + r_c}
$$

where:
- $r_a$: aerodynamic resistance
- $r_b$: quasi-laminar boundary layer resistance
- $r_c$: canopy resistance

### Canopy Resistance

$$
\frac{1}{r_c} = \frac{1}{r_s + r_m} + \frac{1}{r_{lu}} + \frac{1}{r_{dc}} + \frac{1}{r_{cl}}
$$

where:
- $r_s$: stomatal resistance
- $r_m$: mesophyll resistance
- $r_{lu}$: resistance of upper canopy
- $r_{dc}$: resistance of lower canopy/ground
- $r_{cl}$: resistance of exposed surfaces

## Parameter Tables

### Vegetation-Specific Parameters

| Parameter | Conifer | Deciduous | Grass | Units |
|-----------|---------|-----------|--------|-------|
| $V_{c,max}$ (25°C) | 60 | 80 | 40 | μmol/m²/s |
| $J_{max}$ (25°C) | 120 | 160 | 80 | μmol/m²/s |
| $g_1$ | 3.0 | 4.0 | 5.0 | kPa^0.5 |
| $\epsilon_{iso}$ | 0.1 | 10.0 | 0.0 | μg/g/h |
| $c_d$ | 0.15 | 0.20 | 0.10 | - |

### Temperature Response Parameters

| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| $Q_{10,V}$ | 2.0 | - | $V_{c,max}$ temperature sensitivity |
| $Q_{10,J}$ | 1.9 | - | $J_{max}$ temperature sensitivity |
| $Q_{10,R}$ | 2.0 | - | Respiration temperature sensitivity |
| $H_a$ | 72000 | J/mol | Activation energy |
| $H_d$ | 200000 | J/mol | Deactivation energy |
| $\Delta S$ | 650 | J/mol/K | Entropy term |

### Sensitivity Analysis

Key sensitive parameters identified through sensitivity analysis:

1. **Leaf area index** (±20% → ±15% flux change)
2. **Drag coefficient** (±50% → ±25% wind change)
3. **Maximum carboxylation rate** (±20% → ±18% photosynthesis change)

## References

Key scientific references for parameterizations:

1. **Farquhar et al. (1980)**: Photosynthesis model
2. **Ball et al. (1987)**: Stomatal conductance
3. **Guenther et al. (2012)**: Biogenic emissions
4. **Raupach & Thom (1981)**: Canopy turbulence
5. **Dai et al. (2004)**: Two-stream radiation

For implementation details, see the [API Reference](../api/overview.md).
