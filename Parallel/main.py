import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

import geopandas as gpd
import shapely as sp
import vaex as vx
import argparse
import multiprocessing
from datetime import datetime
import h5py
from rich import print, pretty
from h5coro import h5coro, s3driver, filedriver

pretty.install()

from src import *
from time import time as Time

h5coro.config(errorChecking=True, verbose=False, enableAttributes=False)

def dtm():
    return f'[{datetime.now().strftime("%H:%M:%S")}]'

def core_name():
    return int(multiprocessing.current_process().name.split('-')[-1]) % 4

def crs():
    return "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

def process(link, creds, mask, flowgdf, bounds, GPU):
    
    if not link:
        return None
    
    sst = Time()
    
    # list of paths in hdf5
    lasers = []
    beams = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    items = ['land_ice_segments/latitude', 'land_ice_segments/longitude', 'land_ice_segments/h_li', "land_ice_segments/delta_time", "land_ice_segments/fit_statistics/dh_fit_dx", "land_ice_segments/fit_statistics/dh_fit_dy", 'land_ice_segments/ground_track/ref_azimuth','land_ice_segments/atl06_quality_summary']
    path_parts = [[f'/{beam}/{item}' for item in items] for beam in beams]
    names = ['lat', 'lon', 'h_li', 'time', 'dh_fit_dx', 'dh_fit_dy', 'azumith', 'quality', 'rgt']
    
    # join path names
    paths = []
    for prt in path_parts:
        paths.extend(prt)
        
    # read paths 
    h5obj = h5coro.H5Coro(link.split(':')[-1], s3driver.S3Driver, credentials=creds)
    h5obj.readDatasets(paths, block=True)

    # generate coordinate transform
    source_proj4 = '+proj=latlong +datum=WGS84'
    target_proj4 = crs()
    transformer = pyproj.Transformer.from_proj(pyproj.Proj(source_proj4), pyproj.Proj(target_proj4), always_xy=True)
    
    rgt = link.split("_")[2][:4] # extracts the rgt from the filename
    cycle = link.split("_")[2][4:6]
    
    xlim, ylim = bounds
    
    for beam in beams:
        
        lat, lon = h5obj[f'/{beam}/land_ice_segments/latitude'].values, h5obj[f'/{beam}/land_ice_segments/longitude']
        
        # THIS NEEDS TO BE PUT ON A GPU
        x, y = transformer.transform(lon[:], lat[:])
        # find where track escapes bounding box
        idxmin, idxmax = find_nearest(x, xlim[0]), find_nearest(x, xlim[1])
        idymin, idymax = find_nearest(y, ylim[0]), find_nearest(y, ylim[1])
        idxmin, idxmax = min(idxmin, idxmax), max(idxmin, idxmax)
        idymin, idymax = min(idymin, idymax), max(idymin, idymax)
        idmin = max(idxmin, idymin)
        idmax = min(idxmax, idymax)
        
        vdict = {name:h5obj[f"/{beam}/{ext}"].values[idmin:idmax] for name, ext in zip(names, items) if name != 'lat' and name != 'lon'}
        
        vdict['x'] = x[idmin:idmax]
        vdict['y'] = y[idmin:idmax]
        
        # find if the satellite is ascending or descending in latitude
        vdict['ascending'] = calc_direction(lat[idmin:idmax])
        # multiply the along track slope to account for this
        vdict['dh_fit_dx'] *= vdict['ascending']
        
        lasers.append(pd.DataFrame(vdict))
        
    print(f"{dtm()} - CORE: {core_name()} - [bold]Imported[/bold] rgt {rgt}-{beams} in: [bright_cyan]{round(Time()-sst, 4)}[/bright_cyan]s")
    
    # clip lasers by masks, and remove poor quality points
    lst = Time()
    lasers, beams = reduce_dfs(lasers, beams, mask)
    print(f"{dtm()} - CORE: {core_name()} - Clipped rgt: {rgt}-{beams} in [bright_cyan]{round((Time()-lst)*(10**3), 1)}[/bright_cyan]ms")
    lst = Time()
    lasers, beams = dist_azumith(lasers, beams)
    print(f"{dtm()} - CORE: {core_name()} - Calc'd atd & azumith for {rgt}-{beams} in [bright_cyan]{round((Time()-lst)*(10**3), 1)}[/bright_cyan]ms")
    lst = Time()
    lasers, beams = comp_flowslope(lasers, beams, flowgdf)
    print(f"{dtm()} - CORE: {core_name()} - Calc'd flowslope for {rgt}-{beams} in [bright_cyan]{round((Time()-lst), 1)}[/bright_cyan]s")
    
    for laser, beam in zip(lasers, beams):
        fig, ax = plt.subplots(1, 1, figsize = (20, 4))
        ax.plot(laser['along_track_distance'], laser["slope"])
        ax.set_title("dh_fit_dx along track")
        ax.set_ylabel("Flowslope (m/m))")
        ax.set_xlabel("Along track dist (non-datumed) (km)")
        plt.savefig(f"out{core_name()}-{rgt}-{beam}.png")
    
    return {"rgt":rgt, "cycle":cycle, "lasers":lasers}


def find_gline(dct, gline):
    
    try :
        rgt, cycle, lasers = dct['rgt'], dct['cycle'], dct['lasers']
    except TypeError:
        return None
    
    # iterate through each laser, finding intersection of each track and adding offset
    ib = []
    for laser in lasers:
        xs, ys = laser['x'].values, laser['y'].values
        xys = np.vstack((xs, ys)).T
        gdf_points = gpd.GeoDataFrame({"geometry":[Point(x,y) for x,y in zip(xs, ys)]}, crs=crs())
        pi, li = get_intersections([LineString(xys), gline.iloc[0]])
        laser = offset_by_intersect(laser, gdf_points, pi)
        if type(laser) == type(None):
            continue
        laser = interp_clean_single(laser, 5000)
        if type(laser) == type(None):
            continue
        
        # remove track if too many nan
        fill = 0.66
        nonnan = np.count_nonzero(~np.isnan(laser["slope-filt"]))
        if nonnan < len(laser["slope-filt"]) * fill:
            if debug:
                print(f"Nan count of track too high {nonnan} < {len(track['slope-filt']) * fill}")
            return None

        # take deriv
        laser = deriv_on_gpd(laser)

        slope_breaks = [0.003, 0.002, 0.001, 0.0005, 0.0001]
        for thresh in slope_breaks:
            # find peaks
            peak_dists = find_peaks(laser, thresh, 400)
            if len(peak_dists) > 3:
                break

                
        debug=False
                
        # filt peaks
        peakdf = filt_peaks(peak_dists, laser, debug=debug)
        
        # remove certain peak predictions based on individual values (NOT QUALITY SCORE!)
        final_peaks = reduce_peakdf(peakdf)

        # round peak values for display
        round_peaks = [str(round(peak, 2)) for peak in final_peaks["loc"]]
        if debug:
            print(f"REDUCED PICKS:\n {', '.join(round_peaks)}")


        if len(final_peaks) > 0:

            # choose the peak with the highest quality score and output.
            final_peaks = final_peaks.sort_values(by=['qs'], ascending=False)

            peak = final_peaks.iloc[0]["loc"]
            qs = final_peaks.iloc[0]["qs"]
            
            nearest_id = find_nearest(laser["along_dist"], peak)
            ib.append((laser.iloc[nearest_id]["x"], laser.iloc[nearest_id]["y"]))
            print(f"Ib @ {xys[-1]}")

            if debug:
                print(f"CHOICE:\n {peak} <-> {ib}")
        
        else:
            
            return None
            
        """
        fig, ax = plt.subplots(1, 1, figsize = (20, 4))
        ax.plot(laser["along_dist"].values, laser["slope-filt"].values)
        ax.set_title("Flowslope along track")
        ax.set_ylabel("Flowslope (m/m))")
        ax.set_xlabel("Along track dist (non-datumed) (km)")
        plt.savefig(f"out{core_name()}.png")
        """
        return ib


def main():

    allibs = []
    
    # import arguments and parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--Cycles", type=str, help = "list of cycles (comma separated)")
    parser.add_argument("-sp", "--SpatialExtent", type=str, help = "filepath to spaital extent polygon")
    parser.add_argument("-gp", "--GPU", action=argparse.BooleanOptionalAction, help = "GPU Enable/disable GPU")
    args = parser.parse_args()

    GPU = args.GPU

    cycles = [int(x) for x in args.Cycles.split(',')]
    print(f"[bold yellow]CYCLES:[/bold yellow] {cycles}")

    spat_ext = args.SpatialExtent
    print(f"[bold yellow]SPATIAL EXTENT:[/bold yellow] [italic]{spat_ext}[/italic]")
    
    cores = multiprocessing.cpu_count()
    print(f"[bold yellow]CORES:[/bold yellow] {cores}")

    # load spatial extent
    mask, xlim, ylim = spatial_extent(spat_ext)

    # load flowslope data into geodataframe
    flowslp = load_flowslope("../Flow/antarctic_ice_vel_phase_map_v01.nc", mask, xlim, ylim, GPU)
    lst = Time()
    gline = load_gline("../BackgroundData/GroundedIce.gpkg", xlim, ylim)
    print(f"{dtm()} - Loaded ref gline in [bright_cyan]{round((Time()-lst), 1)}[/bright_cyan]s")

    # set up query
    region, links = file_query(spat_ext, cycles)

    # establish s3 credentials
    creds = gen_s3()
    
    batchs = len(links) / cores
    if batchs - int(batchs) != 0: batchs = int(batchs) + 1
    missing = 4 * batchs - len(links)
    for i in range(missing): links.append(None)
    
    linkbatchs = np.reshape(links, (batchs, 4))
    cntbatchs = np.reshape(range(1, len(links)+1), (batchs, 4))
    
    multiprocessing.set_start_method("spawn", force=False)

    print(f"{dtm()} - Starting parallel processes")
    
    for batch in range(batchs):
        st = Time()
        print(f"{dtm()} - -= BATCH {batch} =-")
        links, filecounts = linkbatchs[batch], cntbatchs[batch]
        params = [(l, creds, mask, flowslp, (xlim, ylim), GPU) for l in links]

        with multiprocessing.Pool(cores) as pool:
            dcts = pool.starmap(process, params)
            
        print(f"{dtm()} - [italic bold yellow]Hunting for gline intersection[/italic bold yellow]")
        params = [(dct, gline) for dct in dcts]
        with multiprocessing.Pool(cores) as pool:
            ibs = pool.starmap(find_gline, params)
        ibs_c = [i[0] for i in ibs if i is not None]
        allibs.extend(ibs_c)
            
        fig, ax = plt.subplots(1, 1, figsize = (10, 10))
        range_cnt = 20

        ax.set_facecolor("gainsboro")
        gline.plot(ax=ax, color="black")
        
        for ib in allibs:
            ax.scatter(ib[0], ib[1], color="red", s = 3)

        ax.set_title("Map of track, slope break, and ice plain")
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel("EPSG:3031 x (m)")
        ax.set_ylabel("EPSG:3031 y (m)")

        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])

        plt.savefig("out.png")
        
        print(f"{dtm()} - -= BATCH TIME: {round(Time() - st, 4)} =-")
    
if __name__ == "__main__":
    main()