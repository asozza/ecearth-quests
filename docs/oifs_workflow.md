# 🧭 OIFS TL63-EOCENE setup workflow

This guies tries to put together some info on the procedure developed to set up the OIFS TL63-EOCENE configuration for OIFS. The entire machinery is built in a flexible way and it is based on python3, making use of external calls to eccodes and cdo

> The guides assumes you already have a TL63L31 configuration up and running

> A lot of the procedure is taken from the `modify_oifs.ipynb` notebook which builds on the `EoceneOIFS` class 

This guides provide some clarification the generation of the 4 fundamental files required for OIFS. 

- ICMGGINIT: gaussian reduced surface variables
- ICMGGINIUA: gaussian reduced hybrid level variables
- ICMSHINIT: spectral variables
- ICMCLINIT: gaussian reduced 12-month climatological boundary conditition

## Setup the class

First thing, you need to initialize the class. This requires setting up 3 different folders. 

- `idir`: this define the original TL63 data that you want to start from 
- `odir`: this is the path where you want to create the new directory tree. The assumption is that we use a different `datadir` for modification to boundary conditions since current EC-Earth structure does not support such fine tuning
- `herolddir`: DeepMIP protocol relies on the work by Herold et al, 2014, where all the netcdf files are provided to create orography/bathymetry etc.

```python

from eocene_creator import EoceneOIFS
import xarray as xr

eocene = EoceneOIFS(
    idir="/lus/h2resw01/hpcperm/ccpd/ECE4-DATA",
    odir="/lus/h2resw01/hpcperm/ccpd/ECE4-DATA-EOCEENE",
    herold="/lus/h2resw01/hpcperm/ccpd/EPOCHAL/Herold_etal_2014"
)
```

Then using this class, it is possible to prepare a couple of fundamental xarray fileds that we will use later on, which are vegetation, land sea mask and orography. Please notice, that this implies some operations but mostly interpolation toward the TL63 grid. 

```python
vegetation = xr.open_dataset("eocene.prepare_vegetation_zhang()")
landsea =xr.open_dataset(eocene.prepare_herold(flag="landsea_mask"))
orog = xr.open_dataset(eocene.prepare_herold(flag="orography"))
```

## Generate the ICMGG files

The class allow for a pretty simple call

```python
eocene.create_init(landsea=landsea['landsea_mask'], tvl=vegetation['tvl'], tvh=vegetation['tvh'], cvl=vegetation['cvl'], cvh=vegetation['cvh'], sd_orog=sd_orog)
eocene.create_iniua()
```

For the ICMGGINIT file, we set to zero most of the fields, 
we replace the land-sea mask and we inject:
 the high and low vegetation types (tvh, tvl)

 their fractional coverage (cvh, cvl)

 and the subgrid-scale orography standard deviation (sdfor)

The vegetation types and cover are based on a mapping of the Herold Eocene biome dataset using Zhang et al. (2021).

For the ICMGGINIUA file, the most relevant change is setting the specific humidity to 0 in the entire vertical column.

## Generate the ICMSHINIT

The class allow for a pretty simple call again. The idea is to extract the orography from spectral and replace with the the ones actually provided by Herold et al.

Values of the Herold has to be convert from meter to geopotential height

Finally, given that we do not have a specific knowledge of the initial conditions we opted for applying a massive truncation to the spectral varialbes as vorticity and divergence. This is done by erasing all harmonics but the first, and allow the model to start from a sort of cold start. 


```python
eocene.create_sh(orog=orog['orography'])
```

## Generate the climate files

There is again simple script which operatas setting a constant value for the different albedos and the leaf area index. This can be of course improved. 
It works by splitting the variable separately and by calling the `modify_single_grib` file which assign a new variable to the file.

It then repacks the file and set the proper time axis. 

```python
eocene.create_climate()
```
