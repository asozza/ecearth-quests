"""Class to generate Eocene OIFS data."""
import os
import xarray as xr
import numpy as np
from cdo import Cdo
from utils import modify_single_grib, truncate_grib_file, nullify_grib
from utils import modify_value, replace_value
from utils import extract_grid_info, spectral2gaussian
import shutil
from utils import GRIB2, NC4
cdo = Cdo()

class EoceneOIFS():

    def __init__(self, idir, odir, herold, 
                 resolution="TL63L31",
                 startdate="19900101",):
        """
        Initialize the EoceneOIFS class.
        
        Args:
            indir (str): Input directory for OIFS data.
            outdir (str): Output directory for processed data.
            herold (str): Path to the Herold data.
        """
        
        # main dirs
        self.idir = idir
        self.odir = odir
        self.herold = herold

        if not os.path.exists(self.herold):
            raise FileNotFoundError(f"Herold data not found at {self.herold}")
        if not os.path.exists(self.idir):
            raise FileNotFoundError(f"Input data not found at {self.idir}")
        
        # options
        self.resolution = resolution
        self.startdate = startdate

        # interpolate toward oifs
        kind, spectral, _ = extract_grid_info(self.resolution)
        self.gaussian = spectral2gaussian(spectral, kind)

        # defining directories
        self.idir_init = os.path.join(self.idir, 'oifs', resolution, startdate)
        self.odir_init = os.path.join(self.odir, 'oifs', resolution, startdate)

        self.idir_climate = os.path.join(self.idir, 'oifs', resolution, "climate.v020")
        self.odir_climate = os.path.join(self.odir, 'oifs', resolution, "climate.v020")

        self.idir_amip = os.path.join(self.idir, 'amip-forcing')
        self.odir_amip = os.path.join(self.odir, 'amip-forcing')

        # create directories
        for d in [self.odir_init, self.odir_climate, self.odir_amip]:
            if not os.path.exists(d):
                os.makedirs(d)


    def prepare_herold(self, flag=None):
        """
        Create a new topographic variable (land-sea mask, opensea mask, bathymetry, orography, sd_orography)
        from the topography and sd topography data.

        Args:
            data (xarray.Dataset): Topography data.
            flag (str): Type of variable to create. One of:
                        "landsea_mask", "mask_opensea", "bathymetry", "orography", "sd_orography".
        """

        if not flag:
            raise ValueError("Flag must be specified. Options are: landsea_mask, mask_opensea, bathymetry, orography, sd_orography")

        # Paths for topography and sd_topography data
        herold_topo = os.path.join(self.herold, 'herold_etal_eocene_topo_1x1.nc')
        herold_sd_topo = os.path.join(self.herold, 'herold_etal_stddev_subgrid_etopo1_to_eocene_1x1.nc')

        if flag == "sd_orography":
            if not os.path.exists(herold_sd_topo):
                 raise FileNotFoundError(f"Herold data not found at {self.herold}")
            ds = xr.open_dataset(herold_sd_topo)
            ds = ds[["paleo_stddev_subgrid_topo"]].rename({"paleo_stddev_subgrid_topo": "sd_orography"})
            # Fill NaNs with 0 — safest for GRIB encoding
            ds["sd_orography"] = ds["sd_orography"].fillna(0)
            # Add missing attributes
            ds["sd_orography"].attrs["units"] = "m"
            ds["sd_orography"].attrs["long_name"] = "Standard deviation of subgrid-scale paleotopography"
            
        
        else:
             
            # Load datasets
            if not os.path.exists(herold_topo):
                raise FileNotFoundError(f"Herold data not found at {self.herold}")
        
            ds = xr.open_dataset(herold_topo)

            if flag == "landsea_mask":
                 ds["landsea_mask"] = (ds["topo"] > 0).astype(int)

            elif flag == "mask_opensea":
                 ds["mask_opensea"] = (ds["topo"] < 0).astype(int)

            elif flag == "bathymetry":
                 ds["bathymetry"] = -ds["topo"].where(ds["topo"] < 0, 0)

            elif flag == "orography":
                 ds["orography"] = ds["topo"].where(ds["topo"] > 0, 0)

            else:
                 raise ValueError(f"Unknown flag: {flag}")

        filename = os.path.join(self.herold, f"{flag}.nc")

        # Delete 'topo' from the dataset
        if "topo" in ds:
            del ds["topo"]
        
        if os.path.exists(filename):
            os.remove(filename)
        
        # Save to NetCDF
        ds.to_netcdf(filename)

        #Remap where needed
        if flag in ["landsea_mask", "orography", "sd_orography"]:
            if os.path.exists(os.path.join(self.herold, f"{flag}_remap.nc")):
                os.remove(os.path.join(self.herold, f"{flag}_remap.nc"))
            cdo.remapcon(
                f"N{self.gaussian}",
                input=filename,
                output=os.path.join(self.herold, f"{flag}_remap.nc"))
            return os.path.join(self.herold, f"{flag}_remap.nc")
        else:
            return filename

    def create_sic(self, value=0.0):
        """
        Create sea ice data for the Eocene OIFS.
        
        Args:
            value (float): Value to set for sea ice data.
        """
        # Create sea ice data
        icefield = xr.open_dataset(
            os.path.join(self.idir_amip, 'siconcbcs_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc'))
        icefield['siconcbcs'] = icefield['siconcbcs']*value
        outfile = os.path.join(
            self.odir_amip, 'siconcbcs_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc')
        if os.path.exists(outfile):
            os.remove(outfile)
        icefield.to_netcdf(outfile)

    def create_sst(self, A=25, OFFSET=20, T0=5):

        """
        Generate a seasonal SST latitudinal pattern and save it to a netCDF file.
        Broadcasts the pattern to match the original SST data shape.
        The pattern is a cosine function of latitude with a phase offset for seasonal cycle
        Parameters:
        - A: Amplitude of the latitudinal SST pattern.
        - OFFSET: Phase offset in degree for the seasonal pattern.
        - T0: Mean temperature in Celcius.
        """
        inputfile = os.path.join(
            self.idir_amip, 
            'tosbcs_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc'
        )
        outputfile = os.path.join(
            self.odir_amip, 
            'tosbcs_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc'
        )
        sstfield = xr.open_dataset(inputfile)
        sstfield['tosbcs'].shape
        lons = sstfield['lon'].values
        lats = sstfield['lat'].values
        lon2d, lat2d = np.meshgrid(lons, lats)

        # Sinusoidal parameters
        A = 25        # Amplitude in degrees Celsius
        #B = 5      # Amplitude in degrees Celsius
        #beta = 45    # Phase shift in degrees
        k_lat = np.pi / 180    # frequency in lat direction
        #k_lon = np.pi / 180 * 4    # frequency in lon direction
        T0 = 5       # Mean temperature
        OFFSET = 20

        seasonal = np.cos(np.linspace(0,2*np.pi,num=13))[:-1]* OFFSET

        # Create sinusoidal SST pattern
        sst_pattern = []
        for phasing in seasonal:
            sst_pattern.append(A * np.pow(np.cos(k_lat * lat2d + np.pi/180* phasing), 2) + T0) #+ B * np.sin(k_lon * lon2d + np.pi/180*beta) + T0
        sst_stack = np.stack(sst_pattern, axis=0)
        stacksize = sstfield['tosbcs'].shape[0]
        sst_broadcast = np.tile(sst_stack, ((stacksize+11)//12, 1, 1))[:stacksize]
        sstfield['tosbcs'].data = sst_broadcast
        if os.path.exists(outputfile):
            os.remove(outputfile)
        sstfield.to_netcdf(outputfile)

    def create_climate(self):
        """
        Create the ICMCL data for the Eocene OIFS.
        Set the albedo and the LAI to constant values.
        """

        # dictionary for values for each variable to modify
        #match_dict = {
        #    "al": 0.15,
        #    "aluvp": 0.06,
        #    "aluvd": 0.06,
        #    "alnip": 0.06,
        #    "alnid": 0.06,
        #    "lai_lv": 0.,
        #    "lai_hv": 0.,
        #}

        inputfile= os.path.join(self.idir_climate, 'ICMCLECE4')
        outputfile = os.path.join(self.odir_climate, 'ICMCLECE4')   
        

        # split variables to operate on them individually
        #cdo.splitname(
        #    input=inputfile, 
        #    output=os.path.join(self.odir_climate,'ICMCLECE4_temp_'), 
        #    options=GRIB2)

        # use modify grib to set them to the new value
        #for var, new_value in match_dict.items():
        #    modify_single_grib(
        #        inputfile=os.path.join(self.odir_climate, f'ICMCLECE4_temp_{var}.grb'),
        #        outputfile=os.path.join(self.odir_climate, f'ICMCLECE4_mod_{var}.grb'),
        #        variables=[var],
        #        spectral=False,
        #        myfunction=modify_value,
        #        newvalue=new_value
        #    )
        #    os.remove(os.path.join(self.odir_climate, f'ICMCLECE4_temp_{var}.grb'))

        # merge them them together using the order of the filenames
        #variables = list(match_dict.keys())
        #paths = [os.path.join(self.odir_climate, f'ICMCLECE4_mod_{var}.grb') for var in variables]
        #if os.path.exists(os.path.join(self.odir_climate, 'ICMCLECE4_almost')):
        #    os.remove(os.path.join(self.odir_climate, 'ICMCLECE4_almost'))
        #cdo.mergetime(options="-L", input=paths, 
        #            output=os.path.join(self.odir_climate, 'ICMCLECE4_almost'))

        # for some strange reason CDO mess up the time axis. Set it back an absolute time axis
        # to guarantee that files are read in the correct order
        #cdo.settaxis("9999-01-15,00:00:00,1month", input=os.path.join(self.odir_climate, 'ICMCLECE4_almost'),
        #            output=outputfile, options="-a")
        #for path in paths:
        #    os.remove(path)

    def create_sh(self, orog):
        """
        Create the ICMSH data for the Eocene OIFS.
        Truncate to first harmonics all the spectral fields
        Replace the orography with the one from the Herold data.

        Args:
            orog (xarray.DataArray): Orography data to be used for the ICMSH data.
            sd_orog (xarray.DataArray, optional): Standard deviation of orography.
        """
        

        input_spectral = os.path.join(self.idir_init, 'ICMSHECE4INIT')
        output_spectral = os.path.join(self.odir_init, 'ICMSHECE4INIT')
        input_surface = os.path.join(self.idir_init, 'ICMGGECE4INIT')
         

        # erase all orography
        modify_single_grib(
            inputfile=input_spectral,
            outputfile=output_spectral,
            variables='z',
            spectral=True,
            myfunction=replace_value,
            newfield=orog*9.81 #converted to geopotential
        )

        ## ⬅️ Add sd_orography if provided
        #if sd_orog is not None:
        #    modify_single_grib(
        #        inputfile=output_spectral,
        #        outputfile=output_spectral,
        #        variables='sdor',
        #        spectral=False,  # likely gridpoint
        #        myfunction=replace_value,
        #        newfield=sd_orog
        #    )

        # truncate spectral variables to first harmonic (mean value)
        truncate_grib_file(
            inputfile=output_spectral,
            variables=['t','d','vo','lnsp'],
            outputfile=output_spectral,
        )

    def create_init(self, landsea, tvh, tvl, cvh, cvl, sd_orog):
        """
        Create the ICMGGECE4INIT data for the Eocene OIFS.
        Replace landsea mask
        Set subgrid orography to 0, soild type to 3, and vegetation to 0.

        Args:
            landsea (xarray.DataArray): Land-sea mask data to be used for the ICMGE data.
            tvh (xarray.DataArray): Vegetation type data for high vegetation.
            tvl (xarray.DataArray): Vegetation type data for low vegetation.
            cvh (xarray.DataArray): Vegetation cover data for high vegetation.
            cvl (xarray.DataArray): Vegetation cover data for low vegetation.
            sd_orog (xarray.DataArray): Standard deviation of subgrid-scale orography.
        """

        input_surface = os.path.join(self.idir_init, 'ICMGGECE4INIT')
        output_surface = os.path.join(self.odir_init, 'ICMGGECE4INIT')

         # Start by copying the base surface file
        shutil.copy(input_surface, output_surface)


        # Inject sd_orography (sdfor)
        modify_single_grib(
            inputfile=output_surface,
            outputfile=output_surface,
            variables=['sdfor'],
            spectral=False,
            myfunction=replace_value,
            newfield=sd_orog
        )
        # Zero out other subgrid orographic fields
        modify_single_grib(
            inputfile=input_surface,
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['anor', 'isor', 'slor', 'cl', 'chnk'],
            spectral=False,
            myfunction=modify_value,
            newvalue=0.  
        )

        nullify_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['sd']
        )

        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['tvh'],
            spectral=False,
            myfunction=replace_value,
            newfield=tvh
        )

        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['tvl'],
            spectral=False,
            myfunction=replace_value,
            newfield=tvl
        )

        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['cvh'],
            spectral=False,
            myfunction=replace_value,
            newfield=tvh
        )

        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['cvl'],
            spectral=False,
            myfunction=replace_value,
            newfield=tvl
        )

        #erase all subgrid orography
        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            variables=['slt'],
            spectral=False,
            myfunction=modify_value,
            newvalue=1 
        )

        # update the land sea mask
        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT'),
            outputfile=output_surface,
            variables=['lsm'],
            spectral=False,
            myfunction=replace_value,
            newfield=landsea
        )

    def create_iniua(self):
        """
        Create the ICMGGECE4INIUA data for the Eocene OIFS.
        Set the humidity to 0.
        """

        input_levels = os.path.join(self.idir_init, 'ICMGGECE4INIUA')
        output_levels = os.path.join(self.odir_init, 'ICMGGECE4INIUA')

        modify_single_grib(
            inputfile=input_levels,
            outputfile=output_levels,
            variables='q',
            spectral=False,
            myfunction=modify_value,
            newvalue=0.  
        )

    
    def prepare_vegetation(self):
        """"
        Create the ICMGG vegetation data for the Eocene OIFS.
        Replace the vegetation data with the one from the Herold data.
        Set the vegetation type to 0 for all types.
        Perform a mapping from present-day initial conditions
        """


        herold_file = os.path.join(self.herold, "herold_etal_eocene_biome_1x1.nc")
        herold_remap = cdo.remapnn(
            f"N{self.gaussian}", 
            input=herold_file, 
            output=os.path.join(self.herold, "herold_etal_eocene_biome_1x1_N32.nc")
        )
            
        icmgg_file = os.path.join(self.idir_init, "ICMGGECE4INIT")
        if os.path.exists(os.path.join(self.herold, "ICMGG.nc")):
            os.remove(os.path.join(self.herold, "ICMGG.nc"))
        icmgg_remap = cdo.setgridtype(
            "regularnn", 
            input=icmgg_file, 
            output=os.path.join(self.herold ,"ICMGG.nc"),
            options=NC4
        )

        herold = xr.open_dataset(herold_remap)
        icmgg = xr.open_dataset(icmgg_remap)

        biome_dict = {'tvh': {}, 'tvl': {}}
        for vegtype in ["tvh", "tvl"]:
            for i in range(1, 11):
                vegid = icmgg[vegtype].where(herold["prei_biome_hp"] == i).values
                vegid = vegid[~np.isnan(vegid)]
                unique, counts = np.unique(vegid, return_counts=True)
                if unique.size>0:
                    biome_dict[vegtype][i] = int(unique[np.argmax(counts)])
                else:
                    biome_dict[vegtype][i] = None

        eocene_icmgg = icmgg[['tvh', 'tvl', 'cvh', 'cvl']]
        for vegtype in ["tvh", "tvl"]:
            eocene_icmgg[vegtype] = eocene_icmgg[vegtype]*0
            for i in range(1, 11):
                eocene_icmgg[vegtype] = xr.where(
                    herold['eocene_biome_hp'] == i,
                    biome_dict[vegtype][i],
                    eocene_icmgg[vegtype])
        
        for vegtype in ["cvh", "cvl"]:
            eocene_icmgg[vegtype] = eocene_icmgg[vegtype]*0

        if os.path.exists(os.path.join(self.odir_init, "ICMGG_vegetation.nc")):
            os.remove(os.path.join(self.odir_init, "ICMGG_vegetation.nc"))
        eocene_icmgg.to_netcdf(
            os.path.join(self.odir_init, "ICMGG_vegetation.nc")
        )

        return os.path.join(self.odir_init, "ICMGG_vegetation.nc")


    def prepare_vegetation_zhang(self):
        """"
        Alternative method to create the ICMGG vegetation data for the Eocene OIFS.
        Replace the vegetation data with the one from the Herold data.
        Set the vegetation content to 1 for the dominant vegetation type and 0 to the others.
        Perform a mapping using the Zhang et al., 2021 criteria. 
        """


        herold_file = os.path.join(self.herold, "herold_etal_eocene_biome_1x1.nc")
        herold_remap = cdo.remapnn(
            f"N{self.gaussian}", 
            input=herold_file, 
            output=os.path.join(self.herold, "herold_etal_eocene_biome_1x1_N32.nc")
        )

        herold = xr.open_dataset(herold_remap)

        # === Biome to vegetation ID mappings ===
        biome_to_tvh = {
            1: 1, # Tropical forest → Evergreen broadleaf trees
            2: 2, # Warm-temperate forest → Evergreen needleleaf trees
            6: 3, # Temperate forest → Deciduous broadleaf
            7: 4, # Boreal forest → Deciduous needleleaf
        }

        biome_to_tvl = {
            3: 5, # Savanna → Tall grass
            4: 6, # Grassland → Short grass
            5: 7, # Desert → Semidesert
            8: 8, # Tundra → Tundra
            9: 8, # Dry Tundra → Tundra
        }

        # === Create blank data arrays ===
        shape = herold['eocene_biome_hp'].shape
        coords = herold.coords

        tvh = xr.full_like(herold['eocene_biome_hp'], fill_value=0)
        tvl = xr.full_like(herold['eocene_biome_hp'], fill_value=0)
        cvh = xr.zeros_like(tvh)
        cvl = xr.zeros_like(tvl)

        for biome_id in range(1, 10): # assuming biome IDs go from 1 to 9
            mask = herold['eocene_biome_hp'] == biome_id

            if biome_id in biome_to_tvh:
                tvh = xr.where(mask, biome_to_tvh[biome_id], tvh)
                cvh = xr.where(mask, 1.0, cvh)
            elif biome_id in biome_to_tvl:
                tvl = xr.where(mask, biome_to_tvl[biome_id], tvl)
                cvl = xr.where(mask, 1.0, cvl)
            else:
                print(f"Warning: biome {biome_id} not in mapping.")

        # === Assemble final dataset ===
        vegetation_ds = xr.Dataset(
        {
            "tvh": tvh,
            "tvl": tvl,
            "cvh": cvh,
            "cvl": cvl
        },

        )

        # === Save ===
        output_path = os.path.join(self.odir_init, "ICMGG_vegetation.nc")
        if os.path.exists(output_path):
            os.remove(output_path)
        vegetation_ds.to_netcdf(output_path)

        return output_path









