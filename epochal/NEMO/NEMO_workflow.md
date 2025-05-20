# 🧭 NEMO PALEORCA2 setup workflow

This guide outlines the procedure for creating a new NEMO configuration starting from a custom coordinates file, in this specific case the `PALEORCA2` configuration, using the `DOMAINcfg` tool. This short guide assumes you have a coordinates file that can be used as a starting point. In our case, this will the `paleorca` configuration. 

> ⚠️ The workflow assumes NEMO v4.2 and EC-Earth environment.

## Preprocessing

We start from coordinates file provided by the IPSL team, a grid configuration known as `PALEORCA2`. This is coming from NEMO 3.6, which implies it has the old NEMO convention for halo points, and it is made by 182x174 points.

### Adapting the halo

ORCA grids from NEMO 4.2 has two points less in the zonal direction, and one point less in the meridional one. This means that we can "chop" out points to achieve the right number of lon/lat from older initial/boundry conditions file.

ORCA2 was going from 182x149 to 180x148, so PALEORCA2 must go from 182x174 to 180x174.

This can be achieved with cdo

```bash
cdo sethalo,-1,-1 coords_ori.nc coords_halo.nc
```

### Adding the bounds

It is important to add cell bounds to this coordinate file so that we have a reference file for conservative interpolation. This can be achieved with the tool `orca_bounds_new.py` from AQUA available in the specific branch which can create the required cell bounds. The original tool is meant to work with the mesh_mask file, but given that at this stage it is not available yet we modified slightly the code in order to make it work with simply a coordinate file. The adapted version is available at https://github.com/DestinE-Climate-DT/AQUA/blob/orca-coordinates.nc/cli/orca-grids/orca_bounds_new.py 

```bash
./orca_bounds_new.py --stagg F --no-level coords_halo.nc coords_bounds.nc
```

The options for staggering is set to `F` mostly because `paleorca` is an F-grid according to IPSL team, and `--no-level` is set because there is no need to try to produce vertical coordinates (it will not work since there are no such info in the coordinate file) 

### Interpolate the bathymetry

It is then important to produce an interpolated bathymetry to the target coordinate. This is done in our case starting from the eORCA1 bathymetry. Of course, in case of paleoclimate simulations, an appropriate bathymetry file is required. 


>  Please notice that the bathymetry is considered to be positive, so if you use an orography file you need to be flip the signs

>  Another issue is that original eORCA domain file might not have appropriate bounds for remapping. this has to be generated or retrieve for the occasion, the `orca_bounds_new.py` file can be used as well as any other experiment output. 


```bash
cdo remapnn,coords_bounds.nc \
    -setgrid,ORCA1_gridescription.txt \
    -selname,bathy_metry,nav_lon,nav_lat \
    /path/to/eORCA1/domain_cfg.nc \
    PALEORCA_bathy_metry_from_eORCA1.nc
```

## Running the DOMAINcfg

This tool will create the `domain_cfg.nc` file which is fundamental for all the further steps. 

### Compile the tool

Compile the domain tool, available in `sources/nemo-4.2/tools`

```bash
./maketools -m ecearth -n DOMAINcfg
```

### Modify the namelist

There is a `namelist_cfg` namelist which has to be adapted to work with such configuration.

```fortran
&namcfg
   ln_read_cfg = .false.
   cn_fcoord   = "coords_bounds.nc"
   nn_bathy    = 1
   cn_topo     = "PALEORCA_bathy_metry_from_eORCA1.nc"
   cn_bath     = "bathy_metry"
   cn_lon      = "nav_lon"
   cn_lat      = "nav_lat"
/
```

- `ln_read_cfg = .false.` -> tells NEMO to generate a new domain from provided inputs
- `nn_bathy = 1` -> ingest an external bathymetry file (already remapped on the target grid)
- `cn_bath`, `cn_lon`, `cn_lat` -> variable names that must match your bathymetry file

## 🧪 4. Customize Bathymetry (optional)

Modify `src/domain_zgr.F90` to manually open or close specific straits. 

> Please notice that such modification has to be made for both the full and the partial vertical coordinates. The one currently used is the block under `ln_zps`, which is for partial coordinates

> Good policy is to set the if so that it works only for paleorca configuration

Here below, you can find an example modification where Panama and Thailand are connected to the main land and Gibrailtair and Red Sea straits are openend. Of course, this has to be manually verified for each new configuration according to the user needs.

```fortran                                            ! =====================
IF( cp_cfg == "paleorca" .AND. jp_cfg == 2 ) THEN    ! PALEORCA configuration
    !                                                 ! =====================   
    !
    ii0 = 135 ;     ii1 = 135
    ij0 = 131 ;     ij1 = 131
    DO ji = mi0(ii0), mi1(ii1)
      DO jj = mj0(ij0), mj1(ij1)
          bathy(ji,jj) = 284._wp
      END DO
    END DO
    IF(lwp) WRITE(numout,*)
    IF(lwp) WRITE(numout,*) '      paleorca: Gibraltar strait open at i=',ii0,' j=',ij0

    ii0 = 156    ;   ii1 = 156                   ! Bab el mandeb Strait open
    ij0 = 119    ;   ij1 = 119                   ! (Thomson, Ocean Modelling, 1995)
    DO ji = mi0(ii0), mi1(ii1)
      DO jj = mj0(ij0), mj1(ij1)
        bathy(ji,jj) = 137._wp
      END DO
  END DO
  IF(lwp) WRITE(numout,*)
  IF(lwp) WRITE(numout,*) '             paleorca: Bab el Mandeb strait open at i=',ii0,' j=',ij0

  !
  ii0 = 95    ;   ii1 = 97                   ! Panama strait closing
  ij0 = 116    ;   ij1 = 116                 !
  DO ji = mi0(ii0), mi1(ii1)
    DO jj = mj0(ij0), mj1(ij1)
        bathy(ji,jj) = 0._wp
    END DO
  END DO
  IF(lwp) WRITE(numout,*)
  IF(lwp) WRITE(numout,*) '             paleorca: Panama strait closed at i=',ii0,' j=',ij0

                !
  ii0 = 5    ;   ii1 = 5                   ! Thailand closing
  ij0 = 117    ;   ij1 = 117                 !
  DO ji = mi0(ii0), mi1(ii1)
    DO jj = mj0(ij0), mj1(ij1)
        bathy(ji,jj) = 0._wp
    END DO
  END DO
  IF(lwp) WRITE(numout,*)
  IF(lwp) WRITE(numout,*) '             paleorca: Thailand closed at i=',ii0,' j=',ij0
END IF
```

> 🔎 Grid indices in `ncview` differ from Fortran by +2, so adjust accordingly.

Of course, after each modification you need to recompile the tool. 

### Run DOMAINcfg

```bash
mpirun -np 1 ./make_domain_cfg.exe
```

This creates `domain_cfg.nc`. Manual inspection to `bathy_metry` and `top_level` is required to verify that everything worked as expected.


## Generate the `maskutil.nc`

You need to run a short NEMO simulation to generate `mesh_mask.nc`.


### Run a test experiment to produce the `mesh_mask.nc`

There is a test expeirment where nemo is linked and compiled in `sources/nemo-4.2/cfgs/ECEARTH/EXP00`
You will need:

- linking/copying the `domain_cfg.nc` 
- modifying the `namelist_cfg` so that it points to the domain file,
- copy a `namcouple` so that oasis will not crash, a random one should make the job

It is possible then to run a short experiment.

```bash
mpi_run -np 1 ./nemo
```

The simulation will fail but will produce a `mesh_mask.nc` with all the info about the mask

### Create the `maskutil.nc`

Use the provided Python script to generate the final mask:

```bash
python orca2_create.py
```

You can find the script here:
https://github.com/asozza/ecearth-quests/blob/main/epochal/NEMO/orca2_create.py

Please verify that folder structure inside the script fits your needs

## Generate the complementary files

### Runoff

Runoff maps on the right grid can be produce by interpolating as a nearest neighbour the original eORCA1 file. 

```bash
cdo remapnn,$HPCPERM/coords_new.nc \
  -setgrid,/ec/res4/scratch/itas/ece4/alfa/output/nemo/alfa_oce_1m_T_2151-2151.nc \
  runoff-icb_DaiTrenberth_Depoorter_eORCA1_JD.nc \ 
  runoff-icb_DaiTrenberth_Depoorter_paleoORCA2_JD.nc
```

## Weights

There is a tool `sources/nemo-4.2/tools/WEIGHTS` which with 3 different steps allow for for generating the weights for different input files. 

The two fundamental ones are the `Goutorbe_ghflux.nc`, which is the bottom flux, and `zdfiwm_forcing_r720x360.nc` for the internal gravity waves. Also the `woa13-levitus-L31.nc` migth be necessary since it should provide the the right weights for initiliaziation. Of course, for real paleo applications those files must be different and appropriate weights must be generated. 


### Compile the tool

As done for the DOMAINcdg, this can be compiled with b

```bash
./maketools -m ecearth -n WEIGHTS
```

### Edit the namelist

You have to create a namelist, starting from `namelist_bilin`, according to each of the file you need to work with.
Mostly, you will need to modify input and output file, the intermediate steps, and the variables available in the files

```fortran
&grid_inputs
    input_file = '/lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/initial/Goutorbe_ghflux.nc'
    nemo_file = '/lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain/pORCA2pd/domain_cfg.nc'
    datagrid_file = 'remap_gou_grid.nc'
    nemogrid_file = 'remap_porca2_grid.nc'
    method = 'regular'
    input_lon = 'lon'
    input_lat = 'lat'
    nemo_lon = 'glamt'
    nemo_lat = 'gphit'
```

Other options can be edited in lower places of the file, and will improve mostly metadata. 

### Run the tool

You have to run the triplet of commands, always pointing to the same namelist

```bash
./scripgrid.exe namelist_bilin_Goutorbe1
./scrip.exe namelist_bilin_Goutorbe1
./scripshape.exe namelist_bilin_Goutorbe1
```

Then do the same for any other file

## Necessary modules (optional for ATOS HPC2020)

```bash
module reset
module load prgenv/intel intel/2021.4.0 intel-mkl/19.0.5 hpcx-openmpi/2.9.0
module load hdf5-parallel/1.12.2 netcdf4-parallel/4.9.1 ecmwf-toolbox/2023.04.1.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF4_PARALLEL_DIR/lib:$ECCODES_DIR/lib:$HDF5_DIR/lib:$HPCPERM/ecearth4/revisions/main/sources/oasis3-mct-5.2/arch_ecearth/lib
```

> The oasis libraries are required by NEMO and should be modified according to your EC-Erth installation path


## Useful links

- **DOMAINcfg Tool** (for generating the domain):
  https://forge.nemo-ocean.eu/nemo/nemo/-/tree/main/tools/DOMAINcfg
- **NEMO Documentation** (v4.2.0 Appendix):  
  https://zenodo.org/records/6334656
- **EC-Earth GitHub Tools (e.g., bounds script)**:  
  https://github.com/DestinE-Climate-DT/AQUA/blob/orca-coordinates.nc/cli/orca-grids/orca_bounds_new.py
- **Final maskutil generation**:  
  https://github.com/asozza/ecearth-quests/blob/main/epochal/NEMO/orca2_create.py



