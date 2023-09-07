import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

import geopandas as gpd
import shapely as sp
import vaex as vx
import argparse
import multiprocessing
import datetime
import h5py

from src import *
from time import time as Time

def dtm():
    return datetime.datetime.now()

def process(link, s3, GPU, bounds):
    
    if not link:
        return None
    
    st = Time()
    
    # maybe one day try to convert this to vaex
    filedata = h5py.File(s3.open(link,'rb'),'r')
    ancillary_data = filedata["ancillary_data"]
    rgt, region, cycle = ancillary_data["start_rgt"][0], ancillary_data["start_region"][0], ancillary_data["start_cycle"][0]
    lasers = (filedata["gt1l"], filedata["gt1r"], filedata["gt2l"], filedata["gt2r"], filedata["gt3l"], filedata["gt3r"])
    names = ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
    
    # open each laser and load into vaex
    for laser, name in zip(lasers, names):
        make_laser(laser, rgt, region, cycle, name, bounds, GPU)
        
    print(f"time: {Time()-st}")
    # open s3 using vaex
    #df = vx.open(s3.open(link,'rb'))
    

def main():

    # import arguments and parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--Cycles", type=str, help = "list of cycles (comma separated)")
    parser.add_argument("-sp", "--SpatialExtent", type=str, help = "filepath to spaital extent polygon")
    parser.add_argument("-gp", "--GPU", action=argparse.BooleanOptionalAction, help = "GPU Enable/disable GPU")
    args = parser.parse_args()

    GPU = args.GPU

    cycles = [int(x) for x in args.Cycles.split(',')]
    print(f"CYCLES: {cycles}")

    spat_ext = args.SpatialExtent
    print(f"SPATIAL EXTENT: {spat_ext}")
    
    cores = multiprocessing.cpu_count()
    print(f"CORES: {cores}")

    # load spatial extent
    mask, xlim, ylim = spatial_extent(spat_ext)

    # load flowslope data into geodataframe
    flowslp = load_flowslope("../Flow/antarctic_ice_vel_phase_map_v01.nc", mask, xlim, ylim, GPU)

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
        params = [(l, s3, GPU, (xlim, ylim)) for l in links]

        with multiprocessing.Pool(cores) as pool:
            pool.starmap(process, params)
        print(f"{dtm()} - -= BATCH TIME: {Time() - st} =-")
    
if __name__ == "__main__":
    main()