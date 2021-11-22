# https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip5.html
from pyesgf.logon import LogonManager
import os.path
import logging
from pyesgf.search import SearchConnection
import xarray as xr

cordex_domain = sys.argv[1] # CORDEX domain [STR] eg "EUR-44"
openid = sys.argv[2] # your ESGF openid [STR] eg 'https://esgf-data.dkrz.de/esgf-idp/openid/xxxx'
outdir = sys.argv[3] # where to write results [STR] eg /path/to/results


# fixed pameters
project = 'CORDEX' # eg CMIP6
time_frequency = 'day'  # eg 3hr

if not os.path.exists(outdir):
    os.makedirs(outdir)

# to clear logger:
# https://stackoverflow.com/questions/30861524/logging-basicconfig-not-creating-log-file-when-i-run-in-pycharm
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
# change to logging.DEBUG to get debug indo
logging.basicConfig(
    level=logging.INFO,
    filename=outdir + "logfile",
    filemode="w",
    format="%(asctime)-15s %(levelname)-8s %(message)s")
# variables
vars = [
    'hurs',
    'tasmin',
    'tasmax',
    'uas',
    'vas',
    'tas',
    'pr',
    'ps',
    'rsds',
    'rlds']  # ,'sfcWind'
expers = ['rcp26', 'historical', 'rcp85']

# now we do a simple index method to extract bbox (not necessary but makes download smaller)
lonstart = 40
lonend = 60
latstart = 40
latend = 50

# logon manager
lm = LogonManager()
if not lm.is_logged_on():

    lm.logoff()
    lm.is_logged_on()

    lm.logon_with_openid(openid=openid, password=None, bootstrap=True)
    lm.is_logged_on()

    lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
    lm.is_logged_on()


conn = SearchConnection('https://esgf-data.dkrz.de/esg-search', distrib=True)

for exper in expers:
    logging.info("Experimet: " + exper)
    print("Experimet: " + exper)
    for var in vars:
        myvar = var
        logging.info("Variable: " + var)
        print("Variable: " + var)

        ctx = conn.new_context(
            project='CORDEX',
            domain='EUR-44',
            experiment=exper,
            time_frequency=time_frequency,
            variable=var
        )
        #logging.info("Found hits: " + str(ctx.hit_count))
        print("Found hits: " + str(ctx.hit_count))

        # loop through hits
        for hit in range(0, ctx.hit_count):
            #logging.info("retrieving hit: " +str(hit+1)+"/"+str(ctx.hit_count) +" "+var+" "+exper)
            print("retrieving hit: " + str(hit + 1) + "/" +
                  str(ctx.hit_count) + " " + var + " " + exper)
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
                if os.path.isfile(outdir + outname):
                    #logging.info (outname+" already downloaded!")
                    print(outname + " already downloaded!")
                    continue

                try:
                    ds = xr.open_dataset(my_url, decode_times=False)
                    logging.info("Downloaded" + my_url)
                except IOError:
                    #logging.info("Server likely down skipping"+my_url)
                    print("Server likely down skipping" + my_url)
                    continue
                # logging.info(ds)
                #rp = ds[rotate_pole]
                da = ds[myvar]
                #da = da.sel(lat=slice(latS, latN), lon=slice(lonE, lonW))
                da2 = da[:, latstart:latend, lonstart:lonend]

                da2.to_netcdf(outdir + outname)
                logging.info(outname + " done!")
                print(outname + " done!")


