# Data and third-party notices

The MIT licence in [`LICENSE`](LICENSE) applies only to original project code and documentation for which the repository copyright holder has authority to grant that licence. It does not relicense third-party data, derived data, maps, or software dependencies.

The prepared archives in `data/prepared/` and the reference files in `outputs/` are derived data products. Keep this notice and the source attributions when redistributing them.

## Bundled and downloaded data

### BRAN2020 monthly temperature

- Source: Bluelink Ocean Reanalysis – BRAN2020, NCI.
- Citation: Chamberlain, Oke and Fiedler (2021), DOI [10.25914/6009627c7af03](https://doi.org/10.25914/6009627c7af03).
- Terms: the NCI collection metadata identifies the collection as CC BY 4.0. Confirm the current collection terms at the [NCI catalogue record](https://geonetwork.nci.org.au/geonetwork/srv/eng/xml.metadata.get?uuid=f9372_7752_2015_3718) before a public release.
- Derived use here: regional subsetting, interpolation to the common analysis grid, depth selection, and calculation of derived velocity and classification products.

### NOAA OISST v2.1 monthly sea-surface temperature

- Source: NOAA Physical Sciences Laboratory monthly `sst.mon.mean.nc` product.
- Terms: information on NOAA government servers is public domain unless specifically marked otherwise. Do not imply NOAA endorsement or present modified material as official NOAA material.
- Attribution requested by NOAA PSL: “Data provided by the NOAA/OAR/PSL, Boulder, Colorado, USA, from their Web site at https://psl.noaa.gov/”.
- Derived use here: regional subsetting, alignment to the common analysis grid, masking, and use in the surface and prepared products.

### CMIP6 CNRM-CM6-1-HR monthly ocean temperature

- Source: CNRM-CERFACS CNRM-CM6-1-HR, `Omon/thetao`, member `r1i1p1f2`, native `gn` grid; CMIP historical and ScenarioMIP SSP2-4.5 experiments.
- Terms: the [CMIP6 controlled licence record](https://wcrp-cmip.github.io/CMIP6_CVs/docs/CMIP6_source_id_licenses.html) lists CNRM-CM6-1-HR as CC BY 4.0. The [CMIP6 Terms of Use](https://pcmdi.llnl.gov/CMIP6/TermsOfUse/TermsOfUse6-2.html) state that the controlled licence record supersedes older NetCDF-header licence text.
- Citation: [historical](https://doi.org/10.22033/ESGF/CMIP6.4067) and [SSP2-4.5](https://doi.org/10.22033/ESGF/CMIP6.4190).
- Derived use here: regional subsetting, interpolation to the common analysis grid, depth selection, and calculation of derived velocity and classification products.
- Publications using these data must follow the CMIP6 acknowledgement and citation requirements.

### DPIRD marine bioregions

- Source: State of Western Australia, Department of Primary Industries and Regional Development, *DPIRD Bioregions (DPIRD-095)*, 2023.
- Terms: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as recorded in the [WA data catalogue](https://catalogue.data.wa.gov.au/en/dataset/dpird-bioregions).
- Derived use here: the Gascoyne reporting polygon, calculation buffer, and map context.

### Western Australian coastline

- Source: Natural Earth, 1:50m cultural country boundaries.
- Terms: [public domain](https://www.naturalearthdata.com/about/terms-of-use/).
- Derived use here: contextual land geometry in figures.

## Attribution conditions for derived products

CC BY 4.0 material may be shared and adapted provided that appropriate credit, a licence link, and an indication of changes are retained. The prepared archives and published outputs are modified and derived products; the source entries above provide the required attribution and change context.

## Software dependencies

The project installs the version-pinned [VoCC](https://github.com/JorGarMol/VoCC) R package separately at runtime. `renv.lock` records VoCC as AGPL (>= 3); its licence applies to VoCC. This repository does not bundle a copy of the VoCC package source.

Other R packages remain subject to their own licences, which are recorded in `renv.lock` and their respective package distributions. Confirm licence compatibility and organisational approval before any public archival or redistribution.
