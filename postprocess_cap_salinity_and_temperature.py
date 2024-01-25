''' this python script reads the path to a *.hyd file (master file for DFM DWAQ output), and makes 
a mass conservation correction to the associated binary files at the sites of POTW and tributary inflows'''


import sys,os
import numpy as np
import shutil

# get the path to the main hydro input file
input_hyd= os.environ.get('HYDRO_PATH')        # get path to *.hyd file from BASH enviornment variable
run_name = os.path.basename(input_hyd).strip('.hyd')

# get the max and min for salinity and temperature
salt_max= int(os.environ.get('MAXSALT'))
salt_min= int(os.environ.get('MINSALT'))
temp_max= int(os.environ.get('MAXTEMP'))
temp_min= int(os.environ.get('MINTEMP'))

# get the full path and directory name for the hydro inputs 
input_path = os.path.dirname(input_hyd)
input_dir = os.path.basename(input_path)
run_path = os.path.dirname(input_path)

# we will create a new set of dwaq inputs in a directory that has the suffix "_bound_temp_salt"
output_dir = input_dir + '_bound_temp_salt'
output_path = os.path.join(run_path, output_dir)
output_hyd = os.path.join(output_path,run_name + '.hyd')

# make the directory for the output
if not os.path.exists(output_path):
	os.makedirs(output_path)

# read the contents of the input *.hyd file, parse, and return as a dictionary
with open(input_hyd, 'r') as f:
    lines = f.readlines()
hydro_dict = {}
nn = len(lines)
n = 0
while n < nn:
    line = lines[n]
    a = line.split(None,1)
    if len(a)==2:
        key = a[0]
        value = a[1].rstrip('\n').strip("'").strip() # trim trailing whitespace and end carriage
        n = n + 1
    elif len(a)==1:
        key = a[0]
        value = []
        n = n + 1
        endflag = False
        while not endflag:
            line = lines[n]
            if 'end-' + key in line:
                endflag = True
            else:
                value.append(line.rstrip('\n').strip("'").strip())
            n = n + 1
    hydro_dict[key] = value

# compute the number of segments
nseg = int(hydro_dict['number-water-quality-segments-per-layer'])
nlay = int(hydro_dict['number-water-quality-layers'])
npts = nseg*nlay

# copy the hydro file into the output directory, and also copy any other files, except the salinity and temperautre ones 
for key in hydro_dict.keys():
	
	# file names contain "-file" in their key
	if '-file' in key:

		# identify the source and destination for copying
		# the hydrodynamics file is lacking a suffix, needs special treatment
		if (key=='hydrodynamic-file') and (not '.hyd' in hydro_dict[key]):
			src_path = os.path.join(input_path, hydro_dict[key] + '.hyd')
			dst_path = os.path.join(output_path, hydro_dict[key] + '.hyd')
		else:
			src_path = os.path.join(input_path, hydro_dict[key])
			dst_path = os.path.join(output_path, hydro_dict[key])

		# copy everything except temperature and salinity, handle that later
		if not (key=='temperature-file' or key=='salinity-file'):
			print('reading: %s' % src_path)
			print('writing: %s' % dst_path)
			try:
				shutil.copy(src_path, dst_path)
			except:
				print('ERROR copying file, moving on to next one...')

# process the salinity and temperature files
for key in ['temperature-file', 'salinity-file']:

	if key in hydro_dict.keys():

		if len(hydro_dict[key])>0:

			if key=='temperature-file':
				data_max = temp_max
				data_min = temp_min
			elif key=='salinity-file':
				data_max = salt_max
				data_min = salt_min
		
			src_path = os.path.join(input_path, hydro_dict[key])
			dst_path = os.path.join(output_path, hydro_dict[key])
		
			print('reading: %s' % src_path)
			print('writing: %s' % dst_path)
			print('bounding values between %f and %f' % (data_min, data_max))
		
			# flag to catch unexpected end of file
			stopflag = False
			with open(src_path, 'rb') as f_in:
				with open(dst_path, 'wb') as f_out:
					while not stopflag:
						try:
							t = np.fromfile(f_in, dtype = 'i4', count=1)
							data = np.fromfile(f_in, dtype = 'f4', count = npts)
						except:
							stopflag=True
		
						if not stopflag:
		
							# bound temperature by max and min
							data[data<data_min] = data_min
							data[data>data_max] = data_max
		
							# write corrected temperatures to output file
							try:
								print('... writing t_sec=%d' % t[0])
								f_out.write(np.array(t[0], dtype='i4').tobytes())
								f_out.write(data.astype('f4').tobytes())
							except:
								stopflag=True

		else:

			print('%s not found in hydro input file, skipping' % key)

	else:

		print('%s not found in hydro input file, skipping' % key)
