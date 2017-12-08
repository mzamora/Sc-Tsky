# Sc-Tsky
Obtaining sky brightness temperatures for Sc layers based on soundings at NKX

## Stuff needed
For solving the radiative transfer, you'll need to get Streamer and build it: http://stratus.ssec.wisc.edu/streamer/streamer.html

For creating the inputs and interacting with the outputs, you'll need to get the Sc-utils library: https://github.com/mzamora/Sc-utils.

## Input creation
The following programs are used to create input files for the RTM software Streamer:

make_streamer_input_DYCOMS.m creates a profile based on the DYCOMS II RF-01 case. We use it for validating and selecting the options we want to use.

create_all_NKX_inputs.m (uses make_streamer_input.m) Creates inputs taking into account the whole sounding profile and then estimates cloud properties by assuming a well mixed structure

create_all_NKX_inputs_v2.m (uses make_streamer_input_v2.m) Creates inputs taking only the sounding profile above the cloud top (inversion height).

## Streamer files
For running the files in Streamer, I opted for moving them to the testio folder in streamer and running them there. I did it only because F77 can give you trouble if path names are too long.

You also need to put the files_for_streamer/write_usr.f file in the progs folder of Streamer before building the code. Otherwise the post processing functions won't read the output correctly.

So I move them and then I run the batch of all files recursively with something like 
$ for file in ../testio/NKX/*; do ./streamer -s "$file"; done

And then copy them back here (I guess I could have done all that in a script huh). 

## Output reading
For reading the following programs can be used:

process_DYCOMS.m reads DYCOMS different setup outputs and creates plots

read_all_NKX_output.m (uses read_streamer_output.m) Goes through all files and retrieves LW_down at cloud top, which is used to calculate the sky brightness temperature. It also calculates the tropospheric values of q_t (total water mixing ratio), so we can establish the relationship T_sky(tropospheric q_t).

Mónica Zamora Z., 2017. SRAF at UCSD solar.ucsd.edu
