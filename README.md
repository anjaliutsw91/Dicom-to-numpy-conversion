
# Dicom-to-numpy-conversion
Dicom to numpy conversion


## Running

dicom2numpy.py --data_dir ./dicom --numpy_dst_dir ./numpys --search_using_key True

•  numpy_dst_dir : The destination for the numpy files

•  data_dir: Location of dicom files. For each patient, files have to be arranged into CT and RT.

•  search_using_key: True/False depending on if you will be reading the roi_list or would like to search using keywords


<img src="directoryformat.png" align="center" />

## dicom2numpy.py

Required: 

•	"pat_list": The list of patients. This should be the folder name for each patient. Could be MRNs.

•	"alias_list": This is the name you would like to give each patient for the numpy file. This will be the patient identifier. If the names in the patient list suffice, copy the same for alias_list. For eg, MRN 12345 could be 1, the next MRN could be 2. However you would like to name it.

Optional:

• "org_list": The list of roi name( for each roi) Each institution has different naming conventions for Rectum, Bladder and CTV. Please check to make sure correct structures are selected. If not predetermined it can be searched using keywords. The current keywords used in the code can be seen from lines 227-234. Edit these keywords to your convenience.

• "org_alias_list": The list of alias for each roi. It could be the same as org_list of you want to keep original name.
Each institution has different naming conventions for naming OARs as well as tumor volumes. Please check to make sure correct structures are selected. If not predetermined it can be searched using keywords. The current keywords used in the code can be seen from lines 227-234. Edit these keywords to your convenience.
In this code, a structure would be named ROI_X_Y.npy where X is the patient identifier from the alias list and Y is the structure alias for the specific structure


In this code,
Bladder would be named ROI_X_0_5.npy where X is the patient identifier from the alias list.
Rectum would be named ROI_X_0_4.npy where X is the patient identifier from the alias list.
CTV would be named ROI_X_0_7.npy where X is the patient identifier from the alias list.

If you would prefer a different way of saving, it can be changed within the code.

## Citing
If you find this repository useful, please consider citing it using a link to the repo :)
 
