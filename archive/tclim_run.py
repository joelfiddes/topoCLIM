import tmapp_run_evo_FSM as tmapp
import glob

wd = "/home/joel/sim/qmap/GR_data/"
grids = sorted(glob.glob(wd + "sim/*"))

for grid in grids:
    g = grid.split(wd + 'sim/')[1]
    tmapp.main(wd, simdir=g, model='FSM')
