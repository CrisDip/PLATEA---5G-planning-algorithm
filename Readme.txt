How to execute the optimization planning analysis using PLATEA tool.

LIST OF FILES NEEDED:
-PLATEA_main.m
-locate_sensitive_places.m
-candidate_gNBs.m
-micro_gNB.m
-macro_gNB.m
-alpha_heatmaps.m
-kml2struct.m
-ll2utm.m

Inside each file is present a short description about their usage and scope.

Note that the current version of the overall code has been built to execute the optimization planning for the neighborhood of Torrino-Mezzocammino.
For this reason the following files

-pixel_coordinates.mat
-borders_&_sensitive_places.mat
-candidate_sites_coordinates.mat

contain the information related to the area under investigation, from the borders coordinates to the location of sensitive places.
This is to say that if the the user wants to execute the optimization planning using PLATEA in a different area than the current one, the relative files must be updated accordingly.

MINIMUM SYSTEM PREREQUISITES:
Platea algorithm has been originally coded using 
-personal computer equipped with MATLAB R2019b;
-Intel i7 7th generation processor at 2.7GHz;
-8 GB of RAM

To fasten the execution time the whole code has been lately run on
-Dell PowerEdge R230 equipped with Matlab R2019b;
-Intel Xeon E3-1230 v6 3.5GHz processors;
-64 GB of RAM

Note that the Minimum system requisites are enough to make the code correctly run.

As already mentioned, to correctly run the simulator for spefic area the following steps are needed to be executed:

1. Locate the borderes of the area and build a struct containing their coordinates. Save the infpormation inside borders_&_sensitive_places.mat
2. Locate the borders of the sensistive places located inside the considered area and build a struct containing their coordinates. Save the infpormation inside borders_&_sensitive_places.mat
3. Locate the spots that could host a gNB and buil a struct containing their coordinates. Save the information insode candidate_sites_coordinates.mat

(STEPS 1,2,3 can be performed saving a .kml file from Google Maps and then using the function kml2struct to convert it into a struct readable by MATLAB. See the code)

4. Once the overall scenario inspection is completed, update and modify locate_sensitive_places.m and candidate_gNBs.m files accordingly
5. Run PLATEA_main.m (note that the range of the variables alpha_macro and alpha_micro can be modified and adjusted according to the user needs)
6. When PLATEA finishes its analysys, run alpha_heatmaps.m to evaluate which is the best combination in terms of QoS and scenario costs among all the solutions provided by the algorithm


  

