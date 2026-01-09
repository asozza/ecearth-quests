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
from eocene_functions import albedo, compute_slope, vegetation_zhang
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

         # Start by copying the base surface file
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
            inputfile=output_surface,
            outputfile=output_surface,
            variables=['slor'],
            spectral=False,
            myfunction=compute_slope,
            sd_eoc=sd_orog
        )

        # Set anisotropy and soil type to 1
        modify_single_grib(
            inputfile=output_surface,
            outputfile=output_surface,
            variables=['isor', 'slt'],
            spectral=False,
            myfunction=modify_value,
            newvalue=1.  
        )

        # Zero out other subgrid orographic fields
        nullify_grib(
            inputfile=output_surface,
            outputfile=output_surface,
            variables=['sd', 'sdfor', 'anor', 'cl', 'chnk']
        )

        # Modify vegetation variables
        modify_single_grib(
            inputfile=output_surface,
            outputfile=output_surface,
            variables=['tvh', 'tvl', 'cvh', 'cvl'],
            spectral=False,
            myfunction=vegetation_zhang,
            herold_path=self.herold,
            gaussian=self.gaussian
        )

        # update the land sea mask
       # modify_grib_file(
       #     inputfile=current,
       #     outputfile=output_surface,
       #     variables=['lsm'],
       #     spectral=False,
       #     myfunction=replace_value_grib,
       #     newfield=landsea
       # )
        

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

    def aerosols(self):
        """
        Convert Herold Eocene aerosol data (kg/kg) to column-integrated mass (kg/m²)
        on the IFS grid, using the US Standard Atmosphere for density interpolation.
        """

        print("→ Loading Herold aerosol data")
        paleoaerfile = os.path.join(self.herold, "herold_etal_eocene_CAM4_BAM_aerosols.nc")
        aer_paleo = xr.open_dataset(paleoaerfile)

        # Load IFS reference aerosol climatology
        aerfile_ifs_pd = os.path.join(self.idir, 'oifs', 'ifsdata/aerosol_cams_climatology_43R3a.nc')
        aer_ifs_pd = xr.open_dataset(aerfile_ifs_pd)
        aer_ifs_paleo = aer_ifs_pd.copy()

        # Re-bin the Herold aerosol bins
        aer_paleo_newbin = aer_paleo.copy()
        aer_paleo_newbin['DST01'] = aer_paleo['DST01']*0.57
        aer_paleo_newbin['DST02'] = aer_paleo['DST01']*0.39
        aer_paleo_newbin['DST03'] = aer_paleo['DST02'] + aer_paleo['DST03'] + aer_paleo['DST04']
        aer_paleo_newbin['SSLT01'] = aer_paleo['SSLT01']*0.57
        aer_paleo_newbin['SSLT02'] = aer_paleo['SSLT01']*0.39 + aer_paleo['SSLT02'] + aer_paleo['SSLT03']
        aer_paleo_newbin['SSLT03'] = aer_paleo['SSLT04']

        # Check or generate US Standard Atmosphere
        ua_file = os.path.join(self.herold, "us_standard_atmosphere_newlevs.nc")
        if not os.path.exists(ua_file):
            print("→ Generating US Standard Atmosphere data")
            # US Standard Atmosphere data (alt, temp, g, pressure, density, mu)
            data = """
            -1000   21.50   9.810   11.39   1.347   1.821
            0       15.00   9.807   10.13   1.225   1.789
            1000    8.50    9.804   8.988   1.112   1.758
            2000    2.00    9.801   7.950   1.007   1.726
            3000    -4.49   9.797   7.012   0.9093  1.694
            4000    -10.98  9.794   6.166   0.8194  1.661
            5000    -17.47  9.791   5.405   0.7364  1.628
            6000    -23.96  9.788   4.722   0.6601  1.595
            7000    -30.45  9.785   4.111   0.5900  1.561
            8000    -36.94  9.782   3.565   0.5258  1.527
            9000    -43.42  9.779   3.080   0.4671  1.493
            10000   -49.90  9.776   2.650   0.4135  1.458
            15000   -56.50  9.761   1.211   0.1948  1.422
            20000   -56.50  9.745   0.5529  0.08891 1.422
            25000   -51.60  9.730   0.2549  0.04008 1.448
            30000   -46.64  9.715   0.1197  0.01841 1.475
            40000   -22.80  9.684   0.0287  0.003996 1.601
            50000   -2.5    9.654   0.007978 0.001027 1.704
            60000   -26.13  9.624   0.002196 0.0003097 1.584
            70000   -53.57  9.594   0.00052  0.00008283 1.438
            80000   -74.51  9.564   0.00011  0.00001846 1.321
            """
            lines = data.strip().split('\n')
            parsed_data = [list(map(float, line.split())) for line in lines]
            arr = np.array(parsed_data)
            altitude = arr[:,0]
            temperature = 273.15 + arr[:,1]
            g = arr[:,2]
            pressure = 1e2*arr[:,3]
            density = arr[:,4]
            mu = arr[:,5]

            ds = xr.Dataset(
                {
                    'temperature': (['altitude'], temperature),
                    'gravity': (['altitude'], g),
                    'pressure': (['altitude'], pressure),
                    'density': (['altitude'], density),
                    'viscosity': (['altitude'], mu)
                },
                coords={'altitude': (['altitude'], altitude)}
            )

            # log-pressure coordinate & interpolate to Herold levels
            log_pressure = np.log(ds.pressure)
            ds = ds.assign_coords(log_pressure=log_pressure)
            ds_new = ds.swap_dims({'altitude':'log_pressure'})
            new_log_pressure = np.log(aer_paleo.lev)
            ds_interp = ds_new.interp(log_pressure=new_log_pressure, method='linear')
            ds_interp.to_netcdf(ua_file)
            print(f"→ US Standard Atmosphere saved at {ua_file}")
        else:
            ds_interp = xr.open_dataset(ua_file)

        # Map Herold -> IFS variable names
        var_dict = {
            'Sulfates': "SO4",
            'Black_Carbon_hydrophilic': "CB1",
            'Black_Carbon_hydrophobic': "CB2",
            'Mineral_Dust_bin1': "DST01",
            'Mineral_Dust_bin2': "DST02",
            'Mineral_Dust_bin3': "DST03",
            'Organic_Matter_hydrophilic': "OC1",
            'Organic_Matter_hydrophobic': "OC2",
            'Sea_Salt_bin1': "SSLT01",
            'Sea_Salt_bin2': "SSLT02",
            'Sea_Salt_bin3': "SSLT03"
        }

        # Convert to column mass and regrid
        for varname2 in var_dict:
            varname1 = var_dict[varname2]
            print(f"→ Processing {varname1} -> {varname2}")
            new_var = aer_paleo_newbin[varname1] * ds_interp.density
            new_var_rg = regrid_dataset(new_var, regrid_to_reference=aer_ifs_paleo[varname2])
            new_var_rg_int = -new_var_rg.integrate(coord='altitude')
            aer_ifs_paleo[varname2].data = new_var_rg_int.data.astype('float32')

        # Save Eocene aerosol climatology
        output_path = os.path.join(self.odir_init, "aerosol_cams_climatology_43R3a_EOCENO.nc")
        if os.path.exists(output_path):
            os.remove(output_path)
        aer_ifs_paleo.to_netcdf(output_path)
        print(f"→ Eocene aerosol data saved at {output_path}")
        return output_path




    










