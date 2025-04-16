"""Class to generate Eocene OIFS data."""
import os
import xarray as xr
import numpy as np
from cdo import Cdo
from utils import modify_single_grib, truncate_grib_file
from utils import modify_value, replace_value
from utils import extract_grid_info, spectral2gaussian
from utils import GRIB2
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
        Create a new topographic variable (land-sea mask, opensea mask, bathymetry, orography)
        from the topography data.

        Args:
            data (xarray.Dataset): Topography data.
            flag (str): Type of variable to create. One of:
                        "landsea_mask", "mask_opensea", "bathymetry", "orography".
        """

        if not flag:
            raise ValueError("Flag must be specified. Options are: landsea_mask, mask_opensea, bathymetry, orography")

        # Read the topography data
        herold_topo = os.path.join(self.herold, 'herold_etal_eocene_topo_1x1.nc')

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

        # interpolate toward oifs
        kind, spectral, _ = extract_grid_info(self.resolution)
        gaussian = spectral2gaussian(spectral, kind)

        if flag in ["landsea_mask", "orography"]:
            if os.path.exists(os.path.join(self.herold, f"{flag}_remap.nc")):
                os.remove(os.path.join(self.herold, f"{flag}_remap.nc"))
            cdo.remapcon(
                f"N{gaussian}",
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
        match_dict = {
            "al": 0.15,
            "aluvp": 0.06,
            "aluvd": 0.06,
            "alnip": 0.06,
            "alnid": 0.06,
            "lai_lv": 0.,
            "lai_hv": 0.,
        }

        inputfile= os.path.join(self.idir_climate, 'ICMCLECE4')
        outputfile = os.path.join(self.odir_climate, 'ICMCLECE4')   
        
        # split variables to operate on them individually
        cdo.splitname(
            input=inputfile, 
            output=os.path.join(self.odir_climate,'ICMCLECE4_temp_'), 
            options=GRIB2)

        # use modify grib to set them to the new value
        for var, new_value in match_dict.items():
            modify_single_grib(
                inputfile=os.path.join(self.odir_climate, f'ICMCLECE4_temp_{var}.grb'),
                outputfile=os.path.join(self.odir_climate, f'ICMCLECE4_mod_{var}.grb'),
                variables=[var],
                spectral=False,
                myfunction=modify_value,
                newvalue=new_value
            )
            os.remove(os.path.join(self.odir_climate, f'ICMCLECE4_temp_{var}.grb'))

        # merge them them together using the order of the filenames
        variables = list(match_dict.keys())
        paths = [os.path.join(self.odir_climate, f'ICMCLECE4_mod_{var}.grb') for var in variables]
        if os.path.exists(os.path.join(self.odir_climate, 'ICMCLECE4_almost')):
            os.remove(os.path.join(self.odir_climate, 'ICMCLECE4_almost'))
        cdo.mergetime(options="-L", input=paths, 
                    output=os.path.join(self.odir_climate, 'ICMCLECE4_almost'))

        # for some strange reason CDO mess up the time axis. Set it back an absolute time axis
        # to guarantee that files are read in the correct order
        cdo.settaxis("9999-01-15,00:00:00,1month", input=os.path.join(self.odir_climate, 'ICMCLECE4_almost'),
                    output=outputfile, options="-a")
        for path in paths:
            os.remove(path)

    def create_sh(self, orog):
        """
        Create the ICMSH data for the Eocene OIFS.
        Truncate to first harmonics all the spectral fields
        Replace the orography with the one from the Herold data.

        Args:
            orog (xarray.DataArray): Orography data to be used for the ICMSH data.
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
            newfield=orog*9.81 #converted to geopsotential
        )

        # truncate spectral variables to first harmonic (mean value)
        truncate_grib_file(
            inputfile=output_spectral,
            variables=['t','d','vo','lnsp'],
            outputfile=output_spectral,
        )

    def create_init(self, landsea):
        """
        Create the ICMGGECE4INIT data for the Eocene OIFS.
        Replace landsea mask
        Set subgrid orography to 0, soild type to 3, and vegetation to 0.

        Args:
            landsea (xarray.DataArray): Land-sea mask data to be used for the ICMGE data.
        """

        input_surface = os.path.join(self.idir_init, 'ICMGGECE4INIT')
        output_surface = os.path.join(self.odir_init, 'ICMGGECE4INIT')

        # erase all subgrid orography
        modify_single_grib(
            inputfile=input_surface,
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT_temp1'),
            variables=['sdor', 'anor', 'isor', 'slor', 'cl', 'chnk', 'tvh', 'tvl','cvh', 'cvl'],
            spectral=False,
            myfunction=modify_value,
            newvalue=0.  
        )

        # erase all subgrid orography
        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT_temp1'),
            outputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT_temp2'),
            variables=['slt'],
            spectral=False,
            myfunction=modify_value,
            newvalue=3 
        )

        # update the land sea mask
        modify_single_grib(
            inputfile=os.path.join(self.odir_init, 'ICMGGECE4INIT_temp2'),
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









