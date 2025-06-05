#!/usr/bin/python
# -*- coding: utf-8 -*-

# ### Script to compute paleo runoff maps for ece4 starting from orographic slope

import xarray as xr
#from climtools import climtools_lib as ctl
from matplotlib import pyplot as plt
import numpy as np
import xesmf as xe
import argparse

# ### Functions

def follow_point(i, j, flowdir, verbose = False):
    """
    Tracks a point to the next one, following the slope given by flowdir.
    """
    
    compass = 'N NE E SE S SW W NW'.split()

    if verbose: print(flowdir[i,j])

    if flowdir[i,j] == 0:
        if verbose: print('Sea point!')
        return np.nan
    
    dir_ok = compass[int(flowdir[i,j]-1)]
    next = np.array([0, 0])

    if 'N' in dir_ok:
        next += np.array([1, 0])
    if 'S' in dir_ok:
        next += np.array([-1, 0])
    if 'E' in dir_ok:
        next += np.array([0, 1])
    if 'W' in dir_ok:
        next += np.array([0, -1])

    newpo = np.array([i,j]) + next

    newpo_ok = newpo.copy()
    if newpo[0] == 180:
        newpo_ok[0] = 178 # if passing north pole
    if newpo[0] == -1:
        newpo_ok[0] = 1 # if passing south pole
    if newpo[1] == 360:
        newpo_ok[1] = 0 # if passing east boundary
    if newpo[1] == -1:
        newpo_ok[1] = 359 # if passing west boundary
    
    return newpo_ok


def find_sea(i, j, flowdir, verbose = False):
    point = np.array([i,j])
    nupo = point.copy()
    val = flowdir[nupo[0], nupo[1]]

    if val == 0:
        if verbose: print('Already a sea point!')
        return 2, np.nan, np.nan, np.nan

    count = 0
    while val > 0:
        if verbose: print(f'step {count}')
        last = nupo.copy()
        nupo = follow_point(nupo[0], nupo[1], flowdir = flowdir)
        val = flowdir[nupo[0], nupo[1]]
        count += 1

        if count > 100: break
    
    if count > 100:
        if verbose: print('Ended in a loophole!')
        return 1, np.nan, np.nan, np.nan
    else:
        if verbose: print(f'Found sea at ({nupo[0]}, {nupo[1]})')
        return 0, nupo, last, int(flowdir[last[0], last[1]])
    

def track_rivers(flowdir):
    """
    Finds basins and river exit points from a global slope map.
    """
    rivers = []
    rivers_dir = []
    rnf_map = np.zeros((180, 360))-1

    for j in range(360):
        for i in range(180):
            res, po, pocoast, podir = find_sea(i, j, flowdir = flowdir)

            if res == 0:
                # if len(rivers) == 0:
                #     rivers.append(po)
                #     continue
                # if not np.any(np.array(rivers) == po):
                    # rivers.append(po)
                if tuple(po) not in rivers:
                    rivers.append(tuple(po))
                    rivers_dir.append(podir)
                
                rnf_map[i,j] = rivers.index(tuple(po))
                #rnf_map[i,j] = np.argmin(rivers == po)
                
            elif res == 1:
                rnf_map[i,j] = -1
            elif res == 2:
                rnf_map[i,j] = -2
        
    print(f'Found {len(rivers)} rivers!')
    
    return rnf_map, rivers, rivers_dir


def track_rivers_vct(flowdir):
    """
    With arrays, but slower. Don't use!!
    """
    rivers = []
    rivers_dir = []
    rnf_map = np.zeros((180, 360))-1

    for j in range(360):
        for i in range(180):
            res, po, pocoast, podir = find_sea(i, j)

            if res == 0:
                if len(rivers) == 0:
                    rivers.append(po)
                    continue

                if not np.any(np.all(np.array(rivers) == po, axis = 1)):
                    rivers.append(po)
                    rivers_dir.append(podir)
                
                #rnf_map[i,j] = rivers.index(tuple(po))
                rnf_map[i,j] = np.where(np.all(np.array(rivers) == po, axis = 1))[0][0]
                
            elif res == 1:
                rnf_map[i,j] = -1
            elif res == 2:
                rnf_map[i,j] = -2
    
    return rnf_map, rivers, rivers_dir


def plot_basins(rnf_map, rivers = None, not_assigned = None, shuffle = True):
    fig = plt.figure(figsize = (16,12))
    cmap = plt.get_cmap('nipy_spectral').copy()
    cmap.set_under('grey')  # Color for values < vmin

    if shuffle:
        unique_vals = np.unique(rnf_map)
        unique_vals = unique_vals[unique_vals >= 0]
        random_vals = np.random.choice(len(unique_vals), size=len(unique_vals), replace=False)
        randvals = np.select([rnf_map == v for v in unique_vals], random_vals)
        randvals[rnf_map < 0] = -2
        plt.imshow(randvals[::-1, :], cmap=cmap, vmin = 0)
    else:
        plt.imshow(rnf_map[::-1, :], cmap=cmap, vmin = 0)

    if rivers is not None:
        for randv, riv in zip(random_vals, rivers):
            col = cmap(randv/np.max(random_vals))
            plt.scatter(riv[1], 180-riv[0], color = col, s = 20, marker = 'D', edgecolor = 'black')
    
    if not_assigned is not None:
        for po in not_assigned:
            plt.scatter(po[1], 180-po[0], color = 'white', s = 10, marker = 'o', edgecolor = 'black')
    
    return fig


def calc_basins_dim(rnf_map, lat):
    coslat = np.abs(np.cos(np.deg2rad(lat)))
    coslat2d = np.tile(coslat, (360, 1)).T  # Shape

    basins = []
    for ri in range(int(np.max(rnf_map)+1)):
        #basins.append(np.sum(rnf_map == ri))
        basins.append(np.sum(coslat2d[rnf_map == ri]))
        
    basins = np.array(basins)
    return basins


def get_largest_basins(rnf_map, rivers, rivers_dir, lat, basin_thres = 50):
    basins = calc_basins_dim(rnf_map, lat = oroslope.latitude)

    print(f'Selecting {np.sum(basins > basin_thres)} largest basins!')

    big_basins = list(np.where(basins > basin_thres)[0])

    big_rivers = [rivers[i] for i in big_basins]
    big_rivers_dir = [rivers_dir[i] for i in big_basins]

    return big_rivers, big_rivers_dir, big_basins


def track_rivers_and_merge(big_rivers, big_rivers_dir, flowdir, riv_thres = 30, dir_thres = 1, use_expanded = True, weight_for_lat = False, lat = None):
    """
    Tracks rivers and merges small basins to get a reasonable global number of basins.

    riv_thres = 30 # max distance of river to merge
    dir_thres = 1 # max diff in direction
    use_expanded = True # if this is True, it will look not only at the original river estuary, but also at the connected points (slower, but should produce less river basins and ideally more connected)
    weight_for_lat = False # if this is True, the distance is computed weighting for latitude

    """
    coslat = np.abs(np.cos(np.deg2rad(lat)))
    coslat2d = np.tile(coslat, (360, 1)).T  # Shape

    rivers_v2 = [] + big_rivers
    rivers_v2_dir = [] + big_rivers_dir
    rivers_merged = dict()
    for ii, cos in enumerate(rivers_v2):
        rivers_merged[ii] = [np.array(cos)]
    # coastlines = []

    rnf_map_merged = np.zeros((180, 360))-1
    for j in range(360):
        for i in range(180):
            res, po, pocoast, podir = find_sea(i, j, flowdir = flowdir)

            if res == 0:
                if tuple(po) in big_rivers: 
                    print('This is a big river! keeping as is')
                    if tuple(po) not in rivers_v2:
                        rivers_v2.append(tuple(po))
                        rivers_v2_dir.append(podir)
                
                    rnf_map_merged[i,j] = rivers_v2.index(tuple(po))
                else:
                    print('Small river, checking nearby')
                    if tuple(po) not in rivers_v2:
                        ## now for merged
                        # find a river close-by
                        # if the direction is similar (podir +/- 1), then merge
                        # if the dir is different or no river close-by, add new river
                        # rivers merged is a dict of lists

                        # Adding new points assigned to each basin to obtain connected regions
                        if use_expanded:
                            rivers_expanded = rivers_v2.copy()
                            rivers_dir_expanded = rivers_v2_dir.copy()
                            rivers_idx = np.arange(len(rivers_v2))
                            for ri in rivers_merged:
                                rivers_expanded += rivers_merged[ri][1:]
                                rivers_dir_expanded += (len(rivers_merged[ri])-1)*[rivers_v2_dir[ri]]
                                rivers_idx = np.append(rivers_idx, (len(rivers_merged[ri])-1)*[ri])

                            riv_ok = np.abs(np.array(rivers_dir_expanded) - podir) <= dir_thres
                            ok_inds = rivers_idx[riv_ok]

                            sqdist = np.sum((np.array(rivers_expanded)[riv_ok] - np.array(po))**2, axis = 1)
                        else:
                            riv_ok = np.abs(np.array(rivers_v2_dir) - podir) <= dir_thres
                            rivers_idx = np.arange(len(rivers_v2))
                            ok_inds = rivers_idx[riv_ok]

                            sqdist = np.sum((np.array(rivers_v2)[riv_ok] - np.array(po))**2, axis = 1)

                        closest = int(ok_inds[np.argmin(sqdist)])

                        if weight_for_lat:
                            dist = coslat2d[po[0], po[1]] * np.sqrt(np.min(sqdist))
                        else:
                            dist = np.sqrt(np.min(sqdist))

                        print(dist, closest)
                        if dist < riv_thres:
                            print('Found small river nearby!')
                            rnf_map_merged[i,j] = closest
                            if not np.any(rivers_merged[closest] == po):
                                rivers_merged[closest].append(po)
                        else:
                            print('No small river nearby, adding one')
                            # not found a river close-by, adding new one
                            rivers_v2.append(tuple(po))
                            rivers_v2_dir.append(podir)

                            ind = rivers_v2.index(tuple(po))
                            rivers_merged[ind] = [po]
                            rnf_map_merged[i,j] = ind

                    else:
                        # tuple(po) already in rivers_v2
                        rnf_map_merged[i,j] = rivers_v2.index(tuple(po))

            elif res == 1:
                rnf_map_merged[i,j] = -1
            elif res == 2:
                rnf_map_merged[i,j] = -2
        
        print(f'Found {len(rivers_v2)} rivers with {len(big_rivers)} big rivers!')

    return rnf_map_merged, rivers_merged, rivers_v2, rivers_v2_dir


def drop_assign_fill(rnf_map_merged, rivers_v2, rivers_v2_dir, lat = None, basin_thres = 50):
    """
    After the merge, some basins are still below threshold. This routine drops the small basins and assigns the points to one of the closest basins (takes a random neighboring point). The routine goes in loops to assign all points.
    """

    big_rivers_2, big_rivers_dir_2, ok_indx = get_largest_basins(rnf_map_merged, rivers_v2, rivers_v2_dir, lat = lat, basin_thres = basin_thres)

    for ind in range(len(rivers_v2)):
        if ind not in list(ok_indx):
            print(f'{ind}: small basin, deassigning..')
            rnf_map_merged[rnf_map_merged == ind] = -1


    rnf_map_filled = rnf_map_merged.copy()
    nla, nlo = rnf_map_merged.shape

    missing = True
    i_count = 0
    while missing and i_count < 20:
        print(f'We are at {i_count} round of assignments')
        i_count += 1
        missing = False
        not_assigned = []
        for i in range(nla):
            for j in range(nlo):
                if rnf_map_filled[i,j] == -1:
                    print('To be assigned')
                    # Search for neighboring points

                    #print(i-1, i+2, j-1, j+2)
                    ini_i = i-1
                    ini_j = j-1
                    if ini_i < 0: ini_i = 0
                    if ini_j < 0: ini_j = 0
                    neighbors = rnf_map_filled[ini_i:i+2, ini_j:j+2].flatten()
                    print(neighbors)

                    if np.any(neighbors >= 0):
                        # assign closest non-zero
                        new_ind = neighbors[np.where(neighbors >= 0)[0][0]]
                        print(f'Assigning to {new_ind}')
                        rnf_map_filled[i, j] = int(new_ind)
                    else:
                        missing = True
                        not_assigned.append((i,j))

    return rnf_map_merged, not_assigned


def iter_track(flowdir, basin_thres = 50, riv_thres = 30, dir_thres = 1, lat = None, n_iter = 5):
    rnf_map, rivers, rivers_dir = track_rivers(flowdir)

    big_rivers, big_rivers_dir, big_basins = get_largest_basins(rnf_map, rivers, rivers_dir, lat = lat, basin_thres = basin_thres)

    rnf_map_merged, rivers_merged, rivers, rivers_dir = track_rivers_and_merge(big_rivers, big_rivers_dir, flowdir = flowdir, riv_thres = riv_thres, dir_thres = dir_thres, use_expanded = True, weight_for_lat = False, lat = lat)

    ### Now, loop: keep only largest, assign remaining according to a fixed rivers_merged (not updating and not adding new rivers)
    rnf_map_merged_final, not_assigned = drop_assign_fill(rnf_map_merged, rivers, rivers_dir, lat = lat, basin_thres = basin_thres)

    return rnf_map_merged_final, rivers_merged, rivers, not_assigned


def create_basin_data(output_file, ds_target, drainage_basin_id, arrival_point_id, calving_point_id, lon_source, lat_source, method='nearest_s2d'):
    """
    Creates drainage_basin_id, arrival_point_id, and calving_point_id arrays
    based on an existing target dataset and saves to a NetCDF file.
    
    Parameters:
    -----------
    ds_target : xarray.Dataset
        Target dataset with desired grid dimensions (512x256)
    output_file : str
        Path to save the output NetCDF file
    method : str, optional
        Regridding method (default 'nearest_s2d')
    
    Returns:
    --------
    xarray.Dataset
        Dataset containing the three data variables at target resolution
    """
    
    # Create source dataset
    ds_source = xr.Dataset(
        data_vars={
            'drainage_basin_id': (['lat', 'lon'], drainage_basin_id),
            'arrival_point_id': (['lat', 'lon'], arrival_point_id),
            'calving_point_id': (['lat', 'lon'], calving_point_id)
        },
        coords={
            'lon': lon_source,
            'lat': lat_source
        }
    )
    
    # Create regridder
    regridder = xe.Regridder(ds_source, ds_target, method)
    
    # Apply regridding and create target dataset
    result = xr.Dataset(
        data_vars={
            'drainage_basin_id': (['lat', 'lon'], regridder(ds_source.drainage_basin_id).data.astype(np.int32)),
            'arrival_point_id': (['lat', 'lon'], regridder(ds_source.arrival_point_id).data.astype(np.int32)),
            'calving_point_id': (['lat', 'lon'], regridder(ds_source.calving_point_id).data.astype(np.int32))
        },
        coords=ds_target.coords
    )
    
    # Save to NetCDF
    result.to_netcdf(output_file)
    
    return result

# run with: python epochal_ece4_runoff.py -s /abs-path-to-dir/herold_etal_eocene_runoff_1x1.nc -r /abs-path-to-dir/runoff_maps.nc
# output written in runoff_maps_new.nc in the script's folder

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Duplicate job configuration for experiments.")
    parser.add_argument("-s", "--oroslope_file", type=str, required=True, help="Orographic slope file, like herold_etal_eocene_runoff_1x1.nc")
    parser.add_argument("-r", "--runoff_file", type=str, required=True, help="EC-Earth runoff_maps.nc file to update")

    args = parser.parse_args()

    print('Reading oroslope file...')
    oroslope = xr.load_dataset(args.oroslope_file)
    oroslope = oroslope.rename({'xc': 'longitude', 'yc': 'latitude'})
    oroslope = oroslope.assign_coords(longitude = oroslope.longitude[0], latitude = oroslope.latitude[:, 0])
    flowdir = oroslope['RTM_FLOW_DIRECTION'].values

    print('Launching calc!')
    rnf_map_merged_final, rivers_merged_final, rivers_end_point, not_assigned = iter_track(flowdir, lat = oroslope.latitude)

    if len(not_assigned) > 0:
        print('WARNING!! some points have not been assigned:')
        print(not_assigned)

    # Building arrival point id
    arrival_id = np.zeros(rnf_map_merged_final.shape) - 2
    for num in rivers_merged_final:
        for po in rivers_merged_final[num]:
            arrival_id[po[0], po[1]] = num

    ## Setting calving to zero
    calving_id = np.zeros(arrival_id.shape)-2

    rnf_pd = xr.load_dataset(args.runoff_file)
    create_basin_data('runoff_map_new.nc', rnf_pd, rnf_map_merged_final, arrival_id, calving_id, oroslope.longitude, oroslope.latitude)