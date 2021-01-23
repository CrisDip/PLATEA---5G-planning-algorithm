%Name of the file: candidate_gNBs.m
%Script version: v1.0
%Script authors: Luca Chiaraviglio, Cristian Di Paolo

%%This script locates each candidate site to host a
%%gNB in the area under consideration

%Loading the coordinates of each candidate site
load('candidate_sites_coordinates')
k=1;
for i=1:length(total_scenario_struct)
    if total_scenario_struct(i).Description == 'macro'
        macro_struct(i) = total_scenario_struct(i);
    else 
        micro_struct(k) = total_scenario_struct(i);
        k=k+1;
    end
end

candidate_sites_macro = zeros(size(macro_struct,2),2);
candidate_sites_micro = zeros(size(micro_struct,2),2);

for i=1:length(macro_struct)
    macro_struct(i).Name = str2num(macro_struct(i).Name);
    for j=1:size(candidate_sites_macro,2)
        if j == 1
            candidate_sites_macro(i,j) = macro_struct(i).Lat;
        else
            candidate_sites_macro(i,j) = macro_struct(i).Lon;
        end
    end
end

for i=1:length(micro_struct)
    micro_struct(i).Name = str2num(micro_struct(i).Name);
    for j=1:size(candidate_sites_micro,2)
        if j == 1
            candidate_sites_micro(i,j) = micro_struct(i).Lat;
        else
            candidate_sites_micro(i,j) = micro_struct(i).Lon;
        end
    end
end

candidate_sites_macro = ll2utm(candidate_sites_macro);
candidate_sites_micro = ll2utm(candidate_sites_micro);
        