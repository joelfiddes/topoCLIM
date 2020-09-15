# https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip5.html
from pyesgf.logon import LogonManager
import os.path
import logging
from pyesgf.search import SearchConnection
import xarray as xr

#fixed
project='CORDEX'
domain='EUR-44'
time_frequency='day' #3hr
outdir="/home/joel/sim/qmap/raw_cordex/"

print('Downloading' + time_frequency + ' data')
if not os.path.exists(outdir):
    os.makedirs(outdir)

# to clear logger: https://stackoverflow.com/questions/30861524/logging-basicconfig-not-creating-log-file-when-i-run-in-pycharm
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
# change to logging.DEBUG to get debug indo
logging.basicConfig(level=logging.INFO, filename=outdir+"logfile", filemode="w",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
# variables
vars=['hurs', 'tasmin', 'tasmax', 'uas', 'vas', 'tas', 'pr', 'ps', 'rsds','rlds']#,'sfcWind'
expers=['historical','rcp26', 'rcp85']
# define domain here (once figured out rlat/lon translation):
# lonE = 5
# lonW = 11
# latS = 45
# latN = 48

# now we do a simple index method to extract bbox
lonstart=40
lonend=60
latstart=40
latend=50

# logon manager
lm = LogonManager()
if not lm.is_logged_on():

    lm.logoff()
    lm.is_logged_on()

    OPENID = 'https://esgf-data.dkrz.de/esgf-idp/openid/jfiddes'
    lm.logon_with_openid(openid=OPENID, password=None, bootstrap=True)
    lm.is_logged_on()

    lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
    lm.is_logged_on()



conn = SearchConnection('https://esgf-data.dkrz.de/esg-search', distrib=True)

for exper in expers:
    logging.info("Experimet: "+exper)
    print("Experimet: "+exper)
    for var in vars:
        myvar=var
        logging.info("Variable: " +var)
        print("Variable: " +var)

        ctx = conn.new_context(
            project='CORDEX',
            domain='EUR-44',
            experiment=exper,
            time_frequency=time_frequency,
            variable=var
            )
        #logging.info("Found hits: " + str(ctx.hit_count))
        print("Found hits: " + str(ctx.hit_count))

        # f = open("demofile2.txt", "a")
        # f.write("Now the file has more content!")
        # f.close()

        # loop through hits
        for hit in range(0,ctx.hit_count):
            #logging.info("retrieving hit: " +str(hit+1)+"/"+str(ctx.hit_count) +" "+var+" "+exper)
            print("retrieving hit: " +str(hit+1)+"/"+str(ctx.hit_count) +" "+var+" "+exper)
            result = ctx.search()[hit]
            result.dataset_id


            files = result.file_context().search()
            for file in files:

                my_url = file.opendap_url
                if my_url is None:
                    #logging.info("No URL available in hit " + str(hit))
                    print("No URL available in hit " + str(hit))
                    continue

                outname = my_url.split('/')[-1]
                if os.path.isfile(outdir+outname):
                    #logging.info (outname+" already downloaded!")
                    print(outname+" already downloaded!")
                    continue



                try:
                    ds = xr.open_dataset(my_url, decode_times=False)
                    logging.info("Downloaded"+my_url)
                except IOError:
                    #logging.info("Server likely down skipping"+my_url)
                    print("Server likely down skipping"+my_url)
                    continue
                #logging.info(ds)
                #rp = ds[rotate_pole]
                da = ds[myvar]
                #da = da.sel(lat=slice(latS, latN), lon=slice(lonE, lonW))
                da2 = da[:,latstart:latend, lonstart:lonend]

                da2.to_netcdf(outdir+outname)
                logging.info(outname+" done!")
                print(outname+" done!")










# py code to convert between rotataed and normal lonlat from here: TOBE IMPLEMENTED
# https://gis.stackexchange.com/questions/10808/manually-transforming-rotated-lat-lon-to-regular-lat-lon
#
# from math import *
# grid_north_pole_latitude =   39.25
# grid_north_pole_longitude =  -162.0
# north_pole_grid_longitude =  0.0
#
#
# def rotated_grid_transform(lon, lat, option, SP_lon , SP_lat):
#     # lon=9
#     # lat=45
#     # option=2
#     # SP_lon = grid_north_pole_longitude
#     # SP_lat = grid_north_pole_latitude
#     #lon = grid_in[0]
#     #lat = grid_in[1];
#     lon = (lon*pi)/180; # Convert degrees to radians
#     lat = (lat*pi)/180;
#     # SP_lon = SP_coor1;
#     # SP_lat = SP_coor2;
#     theta = 90+SP_lat; # Rotation around y-axis
#     phi = SP_lon; # Rotation around z-axis
#
#     theta = (theta*pi)/180;
#     phi = (phi*pi)/180; # Convert degrees to radians
#
#     x = cos(lon)*cos(lat); # Convert from spherical to cartesian coordinates
#     y = sin(lon)*cos(lat);
#     z = sin(lat);
#
#     if option == 1: # Regular -> Rotated
#
#         x_new = cos(theta)*cos(phi)*x + cos(theta)*sin(phi)*y + sin(theta)*z;
#         y_new = -sin(phi)*x + cos(phi)*y;
#         z_new = -sin(theta)*cos(phi)*x - sin(theta)*sin(phi)*y + cos(theta)*z;
#
#     else:  # Rotated -> Regular
#
#         phi = -phi;
#         theta = -theta;
#
#         x_new = cos(theta)*cos(phi)*x + sin(phi)*y + sin(theta)*cos(phi)*z;
#         y_new = -cos(theta)*sin(phi)*x + cos(phi)*y - sin(theta)*sin(phi)*z;
#         z_new = -sin(theta)*x + cos(theta)*z;
#
#     return(x_new,y_new)
#
# topleftcorner =rotated_grid_transform(lonE,latN,1,grid_north_pole_longitude, grid_north_pole_latitude)
# bottomrightcorner =rotated_grid_transform(lonW,latS,1,grid_north_pole_longitude, grid_north_pole_latitude)
