"""Tests for h5 endpoint."""

import pytest
import h5coro
import webdriver
import earthaccess
import icepyx as ipx
from time import time as Time
import os
import vaex as vx
from rich import print, pretty
pretty.install()

from datetime import datetime
def dtm():
    return f'[{datetime.now().strftime("%H:%M:%S")}]'
    
st = Time()
# set up query object
region = ipx.Query("ATL06", "bungen.gpkg", cycles = [5])
# grab s3 links
gran_ids = region.avail_granules(ids=True, cloud=True)
links = gran_ids[1]
print(f"{dtm()} - Found {len(links)} granules in [bright_cyan]{round(Time() - st, 5)}[/bright_cyan]s")
print(links[0])

auth = earthaccess.login()
credentials = auth.get_s3_credentials(daac="NSIDC")
credentials['aws_access_key_id'] = credentials['accessKeyId']
credentials['aws_secret_access_key'] = credentials['secretAccessKey']
credentials['aws_session_token'] = credentials['sessionToken']

from h5coro import h5coro, s3driver, filedriver

# (2) configure
h5coro.config(errorChecking=True, verbose=False, enableAttributes=False)

# (3) create
st = Time()
beams = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
items = ['land_ice_segments/h_li', "land_ice_segments/delta_time", "land_ice_segments/fit_statistics/dh_fit_dx", "land_ice_segments/fit_statistics/dh_fit_dy", 'land_ice_segments/ground_track/ref_azimuth','land_ice_segments/atl06_quality_summary']
path_parts = [[f'/{beam}/{item}' for item in items] for beam in beams]
names = ['h_li', 'time', 'dh_fit_dx', 'dh_fit_dy', 'azumith', 'quality']
paths = []
for prt in path_parts:
    paths.extend(prt)

for i, link in enumerate(links):
    h5obj = h5coro.H5Coro(link.split(':')[-1], s3driver.S3Driver, credentials=credentials)
    h5obj.readDatasets(paths, block=True)
    for beam in beams:
        #vdict = {name:h5obj[f"/{beam}/{ext}"].values for name, ext in zip(names, items)}
        df = vx.from_dict({name:h5obj[f"/{beam}/{ext}"].values for name, ext in zip(names, items)})
    print(f"{dtm()} - Avg read time: {round((Time()-st)/(i+1), 3)}s", end="     \r")
