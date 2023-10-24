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

pretty.install()

from src import *
from time import time as Time

def dtm():
    return f'[{datetime.now().strftime("%H:%M:%S")}]'

def core_name():
    return int(multiprocessing.current_process().name.split('-')[-1]) % 4

def crs():
    return "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

def process(link, s3, mask, flowgdf, bounds, GPU):
    
    if not link:
        return None
    
    sst = Time()
    
    # maybe one day try to convert this to vaex
    filedata = h5py.File(s3.open(link,'rb'),'r')
    ancillary_data = filedata["ancillary_data"]
    rgt, region, cycle = ancillary_data["start_rgt"][0], ancillary_data["start_region"][0], ancillary_data["start_cycle"][0]
    lasers = (filedata["gt1l"], filedata["gt1r"], filedata["gt2l"], filedata["gt2r"], filedata["gt3l"], filedata["gt3r"])
    names = ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
    
    # open each laser and load into vaex
    df = None
    nlasers = []
    for laser, name in zip(lasers, names):
        lst = Time()
        df = make_laser(laser, rgt, region, cycle, name, bounds, GPU)
        nlasers.append(df)
        st = Time()
    lasers = nlasers
    del nlasers
    print(f"{dtm()} - CORE: {core_name()} - [bold]Imported[/bold] rgt {rgt} in: [bright_cyan]{round(Time()-sst, 4)}[/bright_cyan]s")
    
    # clip lasers by masks, and remove poor quality points
    lst = Time()
    lasers = reduce_dfs(lasers, mask, clpbMsk=False)

    print(f"{dtm()} - CORE: {core_name()} - Clipped rgt: {rgt} in [bright_cyan]{round((Time()-lst)*(10**3), 1)}[/bright_cyan]ms")
    lst = Time()
    lasers = dist_azumith(lasers)
    print(f"{dtm()} - CORE: {core_name()} - Calc'd atd & azumith for {rgt} in [bright_cyan]{round((Time()-lst)*(10**3), 1)}[/bright_cyan]ms")
    lst = Time()
    lasers = comp_flowslope(lasers, flowgdf)
    print(f"{dtm()} - CORE: {core_name()} - Calc'd flowslope for {rgt} in [bright_cyan]{round((Time()-lst), 1)}[/bright_cyan]s")
    
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
        if not laser:
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
    s3 = gen_s3(region)
    
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
        params = [(l, s3, mask, flowslp, (xlim, ylim), GPU) for l in links]

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