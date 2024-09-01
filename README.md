 These scripts set up and run the hydrodynamic solver, DFM, for the SFEI Open Bay model. Much of the setup is automatic! The parts Rusty wrote are totally automatic. The parts Allie wrote require a little more babysitting on the part of the user, but hopefully not too much! 

 The sfb_dfm repo has been migrated to stash and is now found here:
    https://stash.sfei.org/scm/nobm/sfb_dfm
 along with its two submodules:
    https://stash.sfei.org/scm/nobm/sfbay_freshwater
    https://stash.sfei.org/scm/nobm/sfbay_potw
 You can find the "original" sfb_dfm and the original three submodules sfb_dfm_utils, sfbay_freshwater, and sfbay_potw in Rusty Holleman's github account
    https://github.com/rustychris/
 Rusty's version was used to run the original wy2013 simulation. Emma Nuss made some changes, including updating the POTW flows through 2019, and her version is in her github account:
    https://github.com/emmashie
 Note that sfb_dfm_utils has been sucked into the main sfb_dfm repository and is no longer a separate module
 
 You will need access to https://stash.sfei.org/scm/nobm/, and you can email Todd Featherstone (or better yet, make an IT request using Jira) if you need access or have problems.

 Steps for setting up an open bay hydro simulation:

 1) create a folder for the run, e.g. run_folder located at 
        /run_path/run_folder/
    we have been storing our more recent runs in these "run_path" on the new servers:
        /chicagovol1/hpcshared/open_bay/hydro/full_res/
        /boisevol1/hpcshared/open_bay/hydro/full_res/
        /fortcollinsvol1/hpcshared/open_bay/hydro/full_res/
    for example, for a new wy2005 run, on boise, create the following directory:
        /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/
 2) clone sfb_dfm from https://stash.sfei.org/scm/nobm/sfb_dfm into this folder, e.g.
        cd /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/
        git clone https://stash.sfei.org/scm/nobm/sfb_dfm
 3) clone stompy into the same folder:
        cd /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/
        git clone https://github.com/rustychris/stompy
 4) create a folder called "runs" in this same folder, and then create a folder with the name of the run inside that folder, and make a folder called "bc_files" inside that one, e.g., if your run name is "wy2005a" do this:
        cd /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/
        mkdir runs
        cd runs
        mkdir wy2005a
        cd wy2005a 
        mkdir bc_files
    so now you have the following directories in our example:
        /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/
        /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/stompy/
        /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/runs/wy2005a/bc_files/
 5) now navigate inside the sfb_dfm folder and clone the two repositories that 
    contain freshwater inputs from the tributaries and the potw's, respectively. take a moment to make sure that these contain data through the simulation period, e.g. 
        cd /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/
        git clone https://stash.sfei.org/scm/nobm/sfbay_freshwater.git
        git clone https://stash.sfei.org/scm/nobm/sfbay_potw.git 
		
Now you have most of the pieces in place to set up the run. You should take a moment to check that the input files inside all these repositories have data during your intended simulation period. You need to check three things, to start:

 6) check input data to make sure it covers simulation period
	a) check the freshwater inflows file to make sure there are data during the simulation period, e.g.
			/boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/sfbay_freshwater/outputs/sfbay_freshwater.nc
	this data comes primarily from hydrological models at SFEI. for wy2013-wy2017 the data are based on BAHM, and Emma's notes give an overview of that model here 
			https://docs.google.com/document/d/1zcmm4JZ3jDb_MAG8dfY-vgInAC1kAH1N76LW4PV8Q5k/edit#heading=h.emy5l0qu2qlh
	for wy2018-wy2022 we use WDM, and we don't have great documentation for this model yet. talk with Pedro Almodovar, he's the SFEI hydrological modeler, to get more data
	b) check the POTW inflows file to make sure there are data during the simulation period:
			/boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/sfbay_potw/outputs/sfbay_delta_potw.nc
    these data come from a variety of sources, but the vast majority of the data are based on reporting by the POTW's in the annual GAR report that Dave Senn gets from Mike Falk. there's about a one year delay between a given water year and the report availability.

    Warning: make sure the file is actually called "sfbay_delta_potw.nc", and not something
    slightly different like sfbay_delta_potw_Aug2022.nc because this is the file name sfb_dfm will look for
	c) check that the precipitation/evaporation data from CIMIS station 171 in Union City includes the simulation period. the input file is here:
			/boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/sfbay_cimis/union_city-hourly.nc
	if you need to, you can use the scripts in the same folder as this netcdf file to download more data and suck it into the input file. note rusty had this happending automatically but we made it a bit more manual because of major bugginess involved in the automation
    d) to create the salinity and temperature initial condition, and to create a spatially varying secchidepth (for the heat model) we need data from the USGS Peterson Cruise. go to the following website: 
			https://sfbay.wr.usgs.gov/water-quality-database/
    and download data starting from 15 days before the start date of the simulation through the end of the simulation. for a typical single water year simulation, that means downloading data for two calendar years spanning the water year. multiple years of data need to be spliced together by hand. the spliced file should be named wqdata.csv. upload this file to the following folder so sfb_dfm can find it:
			/boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/inputs-static/
	make sure to change the permissions so sfb_dfm can read it:
			cd /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/sfb_dfm/inputs-static/
			chmod ugo+r wqdata.csv
			
The wind and meterological input need to be generated using repositories that live on SFEI's Google Drive. You can either ask Allie to do this (it takes her all of 5 minutes to launch the scripts), or you can try to generate them yourself. To generate them yourself, you either need to use a laptop running Google Drive for Desktop (formerly called Filestream) with the correct folders mounted, or you need to download the whole repository to your computer, and they are big, so you probably don't want to do that. Make sure these inputs are in the UTC time zone like the rest of the model

 7) generate the wind inputs using the following repository:
        1_Nutrient_Share/2_Data_NUTRIENTS/SFEI_Wind/ 
        (link: https://drive.google.com/drive/folders/1e0GGrld8uqjpnBkHCtsQK4-BVl9e1G3G?usp=share_link)
    check the README.txt and Documentation folder to learn how this repo works, or skip ahead and just modify and run the python script /SFEI_Wind/Wind4DFlow-SFB-UTC/generate_amu_amv_4SFB_nearbay_stations_only.py to generate the input files you need. the files are generated in the same folder as this script, and you can pick the name. if you end up with an empty file, you're going to need to download and process more wind data, and you will have to read about how to do that in the documentation folder. note we used to use the file generate_amu_amv_4SFB.py, which used wind from 52 stations around the bay, but in summer 2024 we switched to generate_amu_amv_4SFB_nearbay_stations_only.py which excludes stations that are far from the bay and uses nearest neighbor interpolation instead of linear or natural neighbor, for simplicity. this switch was motivated by improvement of temperature predictions but we did not really properly test whether it improves things, because it was coupled with other changes. nevertheless, it seems like a better idea to only use nearby wind stations and to use a simpler form of interpolation
    
 8) generate the meteorlogical inputs using the following repository:
        1_Nutrient_Share/2_Data_NUTRIENTS/SFEI_Meteo/
        (link: https://drive.google.com/drive/folders/1vph_vQT1CL5BlgomudxDm-58BNfksLQZ?usp=sharing)
    use the script generate_hac_4SFB_curvilinear_station_data_based.py to generate the hac.tem input file (the meteorological forcing). This script replaces the original air temperature (which was from all 52 wind stations) with air temperatures at NDBC stations only (which are on the water). Also replace relative humidity from the gridmet dataset with relative humidity from a set of CIMIS and ASOS stations that are not on the water but are closest to the water. Finally use CIMIS measurements of shortwave radiation to compute cloudiness instead of the daily gridMET dataset, taking care to line up the daily curves in time since sometimes an offset creates crazy cloudiness values. Use simple nearest neighbor interpolation for all three parameters. The data that feeds into this script is all in the SFEI_Wind repository. Note we updated our meteorlogical forcing in summer 2024, and this drastically improves the temperature predictions in our model. The old approach used air temperatures from all 52 wind stations, and used the gridMET dataset to estimate cloudiness and relative humidity. 

  9) upload the wind and meterological forcing files to the following location:
        /run_path/run_folder/runs/run_name/bc_files/
    e.g. in our example, you would upload them into this folder:
        /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/wy2005a/bc_files/
    and make sure to rename them as follows:
        hac.tem = meteorological forcing
        windx.amu = wind forcing, east component
        windy.amv = wind forcing, north component
    This is the location and name where /run_path/run_folder/runs/run_name/FlowFMold_bnd.ext tells DFM they are located. You can take a look at FlowFMold_bnd.ext if you like, it is a text file, and this is where the types of boundary conditions and the paths to their corresponding input files are specified
    
 10) from the linux command line, make sure to change permissions so DFM can read the meteo and wind files, e.g.
        cd /boisevol1/hpcshared/open_bay/hydro/full_res/wy2005/wy2005a/
        chmod -R ugo+r *

 11) if you are running the "old" version of DFM (r52184-opt) you will need to delete the secchidepth parameter from the list of input files, as this parameter is not recognized in the earlier version of the code, and it will crash your simulation. Edit the FlowFMold_bnd.ext file inside the run folder to delete the block of text about secchidepth.

Now we are done with the manual inputs, and we can start the automatic part of the run setup. First, you need to run the python script sfb_dfm.py. This will generate all the remaining input files for your DFlow3D-FM (DFM) run. Next, you need to run DFM on these input files. Finally, if you are planning to use the results of the simulation to run a DWAQ simulation, you are going to need to do two postprocessing steps: first the DWAQ hydro files from the 16 domains in the parallel simulation need to be stitched together, and second, if you are using the old version of the code (r52184-opt), the flow rates coming out of the point sources need to be corrected because this version of DFMleaves them out of the DWAQ hydro files. We have created a series of shell scripts that do all of this for you, but they are very specific to SFEI's servers. If you want to run the model outside SFEI, you will need to find your own way of running sfb_dfm.py, running DFM on your input files, and doing the two postprocessing steps. If you use a newer version of DFM you will not need to correct the point sources, but it doesn't hurt to do that anyway.

Here is how you finally run the model at SFEI: 

There are a series of shell scripts you need to execute in sequence. If you are not planning to use the hydrodynamic model output for DWAQ input you do not need to run step 3 or 4

run_launcher_part0.sh is where you define all the file paths, and the variables defined there are used as input to all the other shell scripts

run_launcher_part1.sh calls sfb_dfm.py to create all the model input files

***Note: if you are using the older version of the code, you will now need to delete the secchidepth input from the .ext file because this version does not accept secchidepth as an input!!! Otherwise your run will crash***

run_launcher_part2.sh partitions the domain and actually runs the DFM solver

run_launcher_part3.sh stitches together the DWAQ input files, which are split across the multidomains used to parallelize the code

run_launcher_part4.sh makes a correction to the DWAQ input files, needed becasue the discharges from the point sources are missing (this error is corrected in newer versions of DFM), and also creates a set of DWAQ input files with temperature and salinity capped at maximum values that are specified in run_launcher_part0.sh

Finally, to create the aggregated grid hydro input for DWAQ, you can run aggregate_hydro.sh 

Optionally, run postprocess_dwaq_binaries_to_netcdf.py to convert the 30 min DWAQ output to something like a map file. You'll need to break the output into short chunks of time, around one month, or server memory will be overloaded and script will crash. Sometimes it crashes anyway and you have to redo it

See Allie's notes about setting up the delft_env anaconda environment here: 
https://docs.google.com/document/d/1M0UWPWKEOPgyxB8YBiivAog91cmQ6KhljN9fR8WRQ2Y/edit#bookmark=id.j7qlzh3zbl0h)

You can run validation scripts that are found here:
	https://stash.sfei.org/scm/nobm/hydro_vaildation_scripts.git

And you may want to aggregate the DWAQ hydro input for use in our fast-running tidally averaged spatially aggregated model. You can find Emma's instructions for doing that here:
    https://docs.google.com/document/d/1KuEs-xHRl-SESOA22cQg1Vkew8vzz5z-ymel3HoBLck/edit#bookmark=id.qnvhwvw92e47
and/or you can borrow and modify the aggregate_hydro.sh script that is included in the sfb_dfm repository to do the job






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

 
