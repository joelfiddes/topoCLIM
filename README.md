# TopoCLIM
Methods to downscale climate timeseries from CORDEX RCM data.

This scheme specifically addresses the need for hillslope scale atmospheric forcing timeseries for modeling the local impact of regional climate change projections on the land surface in complex terrain. The method has a global scope and is able to generate the full suite of model forcing variables required for hydrological and land surface modeling at hourly timesteps. 

It achieves this by utilising the previously published TopoSCALE scheme (Fiddes & Gruber, 2014) to generate a synthetic observation of current climate at hillslope scale, while accounting for a broad range of surface-atmosphere interactions. These synthetic observations are then used to debias (downscale) CORDEX climate variables using the quantile mapping method. A further temporal disaggregation step produces sub-daily fields. This approach has the advantages of other empirical-statistical methods, namely speed of use, while avoiding the need for ground data, which is often limited. It is therefore a suitable method for a wide range of remote regions where ground data is absent, incomplete, or not of sufficient length. The approach is evaluated using a network of high elevation stations across the Swiss Alps and a test application of modelling climate change impacts on Alpine snow cover is given. 

## Getting started

### Installing the `topoCLIM` python module

To begin you'll need a copy of the source code. Either fork the TopoCLIM repository to your own github username, or clone directly.

```{bash}
$ git clone https://github.com/joelfiddes/topoCLIM.git
$ cd topoCLIM
```

It's recommended (but not essential) that you use some sort of python environment manager, such as using the Anaconda distribution and creating an environment (in the code below called "`tclim`"), or using `virtualenv` instead.  This getting started will use Anaconda.

```{bash}
$ conda create -n tclim python
$ conda activate tclim
(tclim)$ pip install -r requirements.txt
```
Install R and R package dependencies:

```{bash}
conda install -c conda-forge r-base
conda install -c conda-forge r-ncdf4
conda install -c omgarcia r-qmap
```

### Running the example
A minimal example for a single point (WFJ) and three climate models is provided here including raw cordex and toposcale downscaled data so it can be run straight out of the box:
```
cd ./topoCLIM/tclim
python run_example.py ../examples/ 1 1
```

### Generating TopoSCALE timeseries

TopoSCALE is a separate python based package available at: https://github.com/joelfiddes/tscaleV2.git

We include a test TopoSCALE timeseries here in "/examples" for running and testing TopoCLIM.

### Downloading CORDEX

To access CORDEX data on the Earth System Grid Federation nodes credentials are required as described here:
https://esgf-data.dkrz.de/user/add/?next=http://esgf-data.dkrz.de/user/add/

CORDEX data is downloaded using the helper script (parameters set in header):
```
python esgf_get.py
```

And postprocessed:
```
python esgf_post.py
```




