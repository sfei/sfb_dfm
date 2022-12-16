 These scripts set up and run the hydrodynamic solver, DFM, for the SFEI Open Bay model
 Much of the setup is automatic! The parts Rusty wrote are totally automatic. The parts Allie
 wrote require a little more babysitting on the part of the user, but hopefully not too much! 

 The sfb_dfm repo has been migrated to stash and is now found here:
    https://stash.sfei.org/scm/nobm/sfb_dfm
 along with its two submodules:
    https://stash.sfei.org/scm/nobm/sfbay_freshwater
    https://stash.sfei.org/scm/nobm/sfbay_potw
 You can find the "original" sfb_dfm and the original three submodules sfb_dfm_utils, 
 sfbay_freshwater, and sfbay_potw in Rusty Holleman's github account
    https://github.com/rustychris/
 Rusty's version was used to run the original wy2013 simulation. Emma Nuss made some
 changes, including updating the POTW flows through 2019, and her version is in her 
 github account:
    https://github.com/emmashie
 Note that sfb_dfm_utils has been sucked into the main sfb_dfm repository and is no 
 longer a separate module
 
 You will need access to https://stash.sfei.org/scm/nobm/, and you can email Todd 
 Featherstone (or better yet, make an IT request using Jira) if you need access or 
 have problems.

 Steps for setting up an open bay hydro simulation:

 1) create a folder for the run, e.g. run_folder located at /run_path/run_folder/
    we have been storing our more recent runs in these "run_path" on the new servers:
        /chicagovol1/hpcshared/open_bay/hydro/full_res/
        /boisevol1/hpcshared/open_bay/hydro/full_res/
        /fortcollinsvol1/hpcshared/open_bay/hydro/full_res/
 2) clone sfb_dfm from https://stash.sfei.org/scm/nobm/sfb_dfm into this folder
        cd /run_path/run_folder/
        git clone https://stash.sfei.org/scm/nobm/sfb_dfm
 3) clone stompy into the same folder:
        git clone https://github.com/rustychris/stompy
 4) now navigate inside the sfb_dfm folder and clone the two repositories that 
    contain freshwater inputs from the tributaries and the potw's, respectively. take 
    a moment to make sure that these contain data through the simulation period
        cd /run_path/run_folder/sfb_dfm/
        clone https://stash.sfei.org/scm/nobm/sfbay_freshwater.git
        clone https://stash.sfei.org/scm/nobm/sfbay_potw.git 
		
Now you have most of the pieces in place to set up the run. You should take a moment to check
that the input files inside all these repositories have data during your intended simulation period!
You need to check three things, to start:

 5) check input data to make sure it covers simulation period
	a) check the freshwater inflows file to make sure there are data during the simulation period:
			/run_path/run_folder/sfb_dfm/sfbay_freshwater/outputs/sfbay_freshwater.nc
	   this data comes primarily from hydrological models at SFEI. for wy2013-wy2017 the data are based
	   on BAHM, and Emma's notes give an overview of that model here 
			https://docs.google.com/document/d/1zcmm4JZ3jDb_MAG8dfY-vgInAC1kAH1N76LW4PV8Q5k/edit#heading=h.emy5l0qu2qlh
	   for wy2018-wy2021 we use WDM, and we don't have great documentation for this model yet. talk with 
	   Tan Zi, he's the modeler, to get more data!
	b) check the POTW inflows file to make sure there are data during the simulation period:
			/run_path/run_folder/sfb_dfm/sfbay_potw/outputs/sfbay_delta_potw.nc
       these data come from a variety of sources, and currently the repository is a bit of a mess. Sienna
	   White worked to clean it up and suck in nutrient input data in addition to flows through wy2019. it is 
	   less messy now, but the part of Rusty's scripts that used to track the data source is broken. the vast 
	   majority of the data are based on reporting by the POTW's in the annual GAR report that Dave Senn gets
	   from Mike Falk. there's about a one year delay between a given water year and the report availability
	c) check that the precipitation/evaporation data from CIMIS station 171 in Union City includes the simulation 
	   period. the input file is here:
			/run_path/run_folder/sfb_dfm/sfbay_cimis/union_city-hourly.nc
	   if you need to, you can use the scripts in the same folder as this netcdf file to download 
	   more data and suck it into the input file. note rusty had this happending automatically but we 
	   made it a bit more manual because of major bugginess involved in the automation
    d) to create the salinity and temperature initial condition, we need data from the USGS Peterson 
       Cruise. go to the following website: 
			https://sfbay.wr.usgs.gov/water-quality-database/
       and download data that spans either side of the start date of the simulation by 15 days. if for some odd 
	   reason you are launching a run near Jan 1, download two years of data and just splice them together by hand. 
	   the file you download will be called wqdata.csv. upload this file to the following folder so sfb_dfm can find it:
			/run_path/run_folder/sfb_dfm/inputs-static/
	   make sure to change the permissions so sfb_dfm can read it:
			cd /run_path/run_folder/sfb_dfm/inputs-static/
			chmod ugo+r wqdata.csv
			
Now it's time to take a break from manual inputs, and start the automatic part of the run setup. You will
be running the sfb_dfm.py script. There is a shell script that sets up all the correct environment variables and 
calls python. Here is how to do it: 

 6) edit run_launcher_part_1.sh inside the sfb_dfm package to specify run name, 
    start time, end time, and flag whether or not you want the boundary conditions
    and sources to be plotted (note plotting boundary conditions also checks for nan's
    so it's a good idea to do this)

 7a) now run run_launcher_part_1.sh to call sfb_dfm.py, which sets up the vast majority of 
	 the model input files. The main input file is the run_name.mdu file, and this file points to 	
     everything else. To execute from command line:
        ./run_launcher_part_1.sh
     you may need to change permissions in order to execute:
        chmod 777 run_launcher_part_1.sh
     This script will set your PYTHONPATH environment variable to point to the copy 
     of stompy installed in the same parent directory as the sfb_dfm package.
	
 7b) if you run into an error running run_launcher_part_1.sh, you may want to run sfb_dfm.py in 
     ipython instead of from a shell script. in this case, you should copy and paste everything 
	 in the run_launcher_part_1.sh script up through right before the following line:
        $PPATH $SFB_DFM_PARENT_PATH/sfb_dfm/sfb_dfm.py
     into your Linux shell. Then if you are on chicago, boise, or fortcollins, activate the anaconda
	 environment we use to run DFM as follows:
        conda activate delft_env
	 (see Allie's notes about setting up the delft_env anaconda environment here: 
	 https://docs.google.com/document/d/1M0UWPWKEOPgyxB8YBiivAog91cmQ6KhljN9fR8WRQ2Y/edit#bookmark=id.j7qlzh3zbl0h)
     launch ipython:
        ipython --pylab
     and run sfb_dfm.py:
        %run sfb_dfm.py

Now we need to take a moment to generate some more input data. The wind and meterological input need to be generated using
repositories that live on SFEI's Google Drive, and you either need to use a laptop running Google Drive for Desktop (formerly
called Filestream) with the correct folders mounted, or you need to download the whole repository to your computer, and they are
big, so you don't want to do that. Make sure these inputs are in the UTC time zone like the rest of the model

 8) generate the wind inputs using the following repository:
		1_Nutrient_Share/2_Data_NUTRIENTS/SFEI_Wind/ 
		(link: https://drive.google.com/drive/folders/1e0GGrld8uqjpnBkHCtsQK4-BVl9e1G3G?usp=share_link)
	check the README.txt and Documentation folder to learn how this repo works, or skip ahead and just 
	modify and run the python script /SFEI_Wind/Wind4DFlow-SFB-UTC/generate_amu_amv_4SFB.py
	to generate the input files you need. the files are generated in the same folder as this script, and you 
    can pick the name. if you end up with an empty file, you're going to need to 
	download and process more wind data, and you will have to read about how to do that in the documentation 
	folder. 
	
 9) generate the meteorlogical inputs using the following repository:
		1_Nutrient_Share/2_Data_NUTRIENTS/SFEI_Meteo/
		(link: https://drive.google.com/drive/folders/1vph_vQT1CL5BlgomudxDm-58BNfksLQZ?usp=sharing)
	sorry, there isn't great documentation for this repo but it's pretty simple. a quick overview is that it 
	contains scripts to download shortwave radiation and relative humidity data from the gridMET dataset, 
	and use those data plus air temperature data from our SFEI_Wind repo, to create an atmospheric forcing input 
	file in the format DFM wants. first look in this folder to make sure you have the data you need: 
	   /SFEI_Meteo/Datasets/gridMET/
	and if you need to download newer data use the metdata_wget.sh shell script (you'll need a terminal emulator
	to do this on a Windows computer, sorry! once you have the data you need, generate your input file using 
	this script:
		/SFEI_Meteo/Meteo4DFlow-SFB-UTC/generate_hac_4SFB_curvilinear.py
	the files are generated in the same folder as this script, and you can pick the name

 10) upload the wind and meterological forcing files to the boundary conditions folder, which 
    should be located here:
        /run_path/run_folder/runs/run_name/bc_files/
    and make sure to rename them as follows:
        hac.tem = meteorological forcing
        windx.amu = wind forcing, east component
        windy.amv = wind forcing, north component
    This is the location and name where /run_path/run_folder/runs/run_name/FlowFMold_bnd.ext tells 
	DFM they are located. You can take a look at FlowFMold_bnd.ext if you like, it is a text file, 
	and this is where the types of boundary conditions and the paths to their corresponding input files are 
	specified
	
 8) from the linux command line, make sure to change permissions so DFM can read the meteo and wind files
        cd /run_path/run_folder/runs/run_name/bc_files/
        chmod ugo+r *
		
Now you have all your input data together and are ready to run the model! the "part 2" run launcher
runs DFM and does some basic postprocessing as well

 9) Edit the user input portion of run_launcher_part_2.sh, and then excecute this script
        ./run_launcher_part_2.sh
    you may need to change permissions to execute:
        chmod 777 run_launcher_part_2.sh

10) Note that if the run executes and generates even partial output, run_launcher_part_2.sh will also 
    automatically perform two postprocessing steps:
        a. Stitch together the multi-domain DWAQ hydro input files into one
        b. Make a correction to the flows in the DWAQ run_name.flo file, necessary because the version of DFM
           we have been using (r52184-opt) leaves out the POTW inflows, violating conservation of mass

Now you are all done, hooray! 

You can run validation scripts that are found here:
	https://stash.sfei.org/scm/nobm/hydro_vaildation_scripts.git

And you may want to aggregate the DWAQ hydro input for use in our fast-running tidally averaged spatially
aggregated model. You can find Emma's instructions for doing that here:
    https://docs.google.com/document/d/1KuEs-xHRl-SESOA22cQg1Vkew8vzz5z-ymel3HoBLck/edit#bookmark=id.qnvhwvw92e47
and/or you can borrow and modify the aggregate_hydro.sh script that is included in the sfb_dfm repository to 
do the job






==================================================
This is what Rusty originally wrote about sfb_dfm:


## What does sfb_dfm.py script do?

* Adds freshwater flows to the model, making sure that the pour points are deep enough (-0.5m) to remain wet at all times.
* Adds POTW flows as point sources at the bed.
* Download tide and Delta outflow data, applying these as additional boundary conditions.
* Update settings in template.mdu to customize paths, time information, output selection.
* Plots boundary conditions (set make_plots = True in sfb_dfm.py)

## Files

**derived**
  Data/files derived from other inputs, but not part of a specific run.
  Currently only used for shapefiles generated from the grid.
  
**inputs-static**
  Files which do not change across runs, are not derived, but are instead
  referenced directly from the MDU or other parts of the model setup.
  
**nudged_features.pli**
  A Delft-style polyline file defining locations of sources.  This is created
  by exporting boundary condition features from Delta Shell.
  
**write_grid_shp.py**
  Short script which writes the shapefile for grid edges, used for loading
  grid representation into GIS.
  
**runs**
  Script-generated simulation setups are in subdirectories below here.
  
**sample_run_dfm**
  Not currently used.  Reference for how to start multiprocessor runs.
  
**sfbay_freshwater**
**sfbay_potw**
  Git submodules holding forcing data for rivers and wastewater discharges.
  
**sfb_dfm.py**
  Main script for generating new runs.
  
**sfei_v19_net.nc**
  Grid.  This grid is modified slightly during the setup of each run based
  on freshwater inputs, and the modified grid is written into the run directory.
  
**template.mdu**
  Template model definition.  Settings which needn't be set dynamically can
  be set here.  Can be tweaked in sfb_dfm.py
  
**update_alviso_bathy.py, plot_sources.py**
  Temporary dev-related scripts for troubleshooting some issues.  Probably
  will be removed down the road.

 
