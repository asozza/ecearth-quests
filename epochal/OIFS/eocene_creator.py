"""Class to generate Eocene OIFS data."""
import os
import xarray as xr
import numpy as np
import xesmf as xe
import shutil
import tempfile
from utils import modify_single_grib, nullify_grib
from utils import modify_value, replace_value, regrid_dataset 
from utils import extract_grid_info, spectral2gaussian
from utils import GRIB2, NC4
from eocene_functions import albedo, compute_slope, vegetation_zhang, aerosols
from cdo import Cdo
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

    def prepare_landsea_mask_present(self, gaussian=48):
        """
        Prepare the present-day land-sea mask from the ICMGGECE4INIT GRIB file.
        Converts to regular Gaussian grid and extracts 'lsm' as an xarray.DataArray.
        """
    
        icmgg_file = os.path.join(self.idir_init, "ICMGGECE4INIT")
    
        # Define temporary NetCDF path
        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
            icmgg_remap = tmp.name

        print(f"→ Reading land-sea mask from {icmgg_file}")
    
        # Convert GRIB to regular Gaussian NetCDF
        cdo.remapnn(
            f"N{self.gaussian}",
            input=f"-setgridtype,regular {icmgg_file}",
            output=icmgg_remap,
            options="-f nc4"
        )
    
        # Open the converted file and extract lsm
        ds = xr.open_dataset(icmgg_remap)
    
        if "lsm" not in ds:
            raise KeyError("Variable 'lsm' not found in the ICMGGECE4INIT file.")
    
        
        # Extract single snapshot of lsm (remove time if exists)
        lsm = ds["lsm"].isel(time=0).squeeze()
        if lsm.dims != ("lat", "lon"):
            lsm = lsm.transpose("lat", "lon")
    
        print(f"Land-sea mask prepared with shape {lsm.shape} and dims {lsm.dims}")
    
        # Clean up temporary file
        os.remove(icmgg_remap)
    
        return lsm

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

    def create_climate(self, lsm_present, landsea):
        """
        Create the ICMCL data for the Eocene OIFS.
        Set the albedo and the LAI to constant values.
        """

        variables = ['al', 'aluvp', 'aluvd', 'alnip', 'alnid', 'lai_lv', 'lai_hv']

        modify_single_grib(
           inputfile=os.path.join(self.idir_climate, "ICMCLECE4"),
           outputfile=os.path.join(self.odir_climate, "ICMCLECE4"),
           variables=variables,
           spectral=False,
           myfunction=albedo,
           lsm_present=lsm_present,
           landsea=landsea  
           ) 


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
        #truncate_grib_file(
        #    inputfile=output_spectral,
        #    variables=['t','d','vo','lnsp'],
        #    outputfile=output_spectral,
        #)

    def create_init(self, landsea, sd_orog):
        """
        Create the ICMGGECE4INIT data for the Eocene OIFS.
        Replace landsea mask
        Modify subgrid orography and vegetation fields, set soil type to 1.

        Args:
            landsea (xarray.DataArray): Land-sea mask data to be used for the ICMGE data.
            sd_orog (xarray.DataArray): Standard deviation of subgrid-scale orography.
        """

        input_surface = os.path.join(self.idir_init, 'ICMGGECE4INIT')
        output_surface = os.path.join(self.odir_init, 'ICMGGECE4INIT')

         # Copying the base surface file
        shutil.copy(input_surface, output_surface)

        # Insert sd_orography (sdor)
        modify_single_grib(
            inputfile=input_surface,
            outputfile=output_surface,
            variables=['sdor'],
            spectral=False,
            myfunction=replace_value,
            newfield=sd_orog
        )

        # Insert slope (slor) computed from sd
        modify_single_grib(
            inputfile=input_surface,
            outputfile=output_surface,
            variables=['slor'],
            spectral=False,
            myfunction=compute_slope,
            sd_eoc=sd_orog
        )

        # Set anisotropy and soil type to 1
        modify_single_grib(
            inputfile=input_surface,
            outputfile=output_surface,
            variables=['isor', 'slt'],
            spectral=False,
            myfunction=modify_value,
            newvalue=1.  
        )

        # Zero out other subgrid orographic fields
        nullify_grib(
            inputfile=input_surface,
            outputfile=output_surface,
            variables=['sd', 'sdfor', 'anor', 'cl', 'chnk']
        )

        # Modify vegetation variables
        modify_single_grib(
            inputfile=input_surface,
            outputfile=output_surface,
            variables=['tvh', 'tvl', 'cvh', 'cvll'],
            spectral=False,
            myfunction=vegetation_zhang,
        )

        # update the land sea mask
        modify_single_grib(
            inputfile=input_surface,
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

    




    










