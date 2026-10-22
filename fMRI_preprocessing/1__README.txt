Purpose: fMRI preprocessing with NORDIC and FIACH
Author: Anna Sadilova
Date: Last updated 21.4.2026

Organisation:
Place scripts into your projects folder - can call on all project data from here
e.g. PROJECTS -> scripts
     PROJECTS -> project -> sub BIDS


Overview of preprocessing scripts:
You can select whether to include NORDIC and FIACH, whether to keep in native space or to convert to MNI, and slice timing
This can be found at the top of the wrapper function:
    run_NORDIC=false
    run_FIACH=true
    MNI_space=false
    slice_timing=true
    
    
------------------------------
BEFORE RUNNING THE FIRST TIME:
------------------------------
This pipeline draws from many different softwares that need to be downloaded - FreeSurfer, FSL, AFNI, Matlab, SPM, ANTS 
     (ants should be automatically downloaded as a python package, so that shouldn't be an issue)



----------------------------------------------------------------------------------------------
** --- Things you'll need to check: --- **
1. That 1_directory.txt has the correct path of project, not scripts
    (IT NEEDS '/' AT THE END e.g. Desktop/Projects/PROJECT/ otherwise the script won't run)

2. Check that all subjects are in a list in 1_subjects.txt

3. If fmap folder is not present, unwarping won't run - fmap needs e1 and e2 phase files!!! 
    Else it will give errors

4. Currently, FIACH is set to output rclean as a non-regressed, non-highpassed dataset
    (only spike interpolated), however as written in AS6, both can be applied via cfg settings
    (For more info, look at AS6)
    


Make sure that all .sh files have been made executable (chmod +x) - this is included in the wrapper function



-------------------------------
To run the preprocessing pipeline, navigate to scripts folder in command line (cd ....) and then call wrapper function by
./wrapper_function_preprocess.sh
--------------------------------



Updated file organisation 21.04.2026:
	Files are now read in fmap and func separately for subfolders 
	- if no subfolders present, it will read data directly and files will be organised as preproc -> all folders here
	- if there are folders in func - e.g. cap1 and cap2, the organisation will be:
		preproc -> T1, checkpoints, cap1 and cap2 folders
		cap1 -> NORDIC,... all other folders
	- folders in preproc are named according to the folders in func
	
	Additionally, checkpoints have been introduced so if there is an error the pipeline, can continue from the same step when issue was fixed 
	(or a step can be skipped altogether if needed - e.g. unwarping fmap folder is empty etc)
	
	
	
	
	
