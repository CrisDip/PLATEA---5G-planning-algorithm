%Name of the file: locate_sensitive_places.m
%Script version: v1.0
%Script authors: Luca Chiaraviglio, Cristian Di Paolo

%%This script has been designed to look for all the
%%sensitive places (schools, churches, public parks) contained in the area
%%under evaluation, with the aim of not placing any candidate site for the 
%%installation of a gNB inside the sensitive areas detected

%Load the struct containing the coordinates of the sensitive places
load('borders_&_sensitive_places.mat')

%setting control variables 
a=1;b=1;c=1;d=1;e=1;f=1;g=1;h=1;l=1;m=1;n=1;o=1;p=1;q=1;

borders = zeros(length(borders_struct),2);
for i=1:length(borders_struct)
    borders_struct(i).Name = str2num(borders_struct(i).Name);
    for j=1:2
        if j == 1
            borders(i,j) = borders_struct(i).Lat;
        else
            borders(i,j) = borders_struct(i).Lon;
        end
    end
end

%shaping torrino map from coordinates
borders = ll2utm(borders);
torrino = polyshape(borders(:,1),borders(:,2));

for i=1:length(park_position_struct)
    park_position_struct(i).Name = str2num(park_position_struct(i).Name);
end

%getting the coordinates of each sensitive place
for i=1:length(park_position_struct)
    if park_position_struct(i).Name == 1
        park1_x(a,1) = park_position_struct(i).Lat;
        park1_y(a,1) = park_position_struct(i).Lon;
        a = a+1;
    elseif park_position_struct(i).Name == 2
        park2_x(b,1) = park_position_struct(i).Lat;
        park2_y(b,1) = park_position_struct(i).Lon;
        b=b+1;
    elseif park_position_struct(i).Name == 3
        park3_x(c,1) = park_position_struct(i).Lat;
        park3_y(c,1) = park_position_struct(i).Lon;
        c=c+1;
    elseif park_position_struct(i).Name == 4
        park4_x(d,1) = park_position_struct(i).Lat;
        park4_y(d,1) = park_position_struct(i).Lon;
        d=d+1;
    elseif park_position_struct(i).Name == 5
        park5_x(e,1) = park_position_struct(i).Lat;
        park5_y(e,1) = park_position_struct(i).Lon;
        e=e+1;
    elseif park_position_struct(i).Name == 6
        park6_x(g,1) = park_position_struct(i).Lat;
        park6_y(g,1) = park_position_struct(i).Lon;
        g=g+1;
    elseif park_position_struct(i).Name == 7
        park7_x(h,1) = park_position_struct(i).Lat;
        park7_y(h,1) = park_position_struct(i).Lon;
        h=h+1;
    elseif park_position_struct(i).Name == 8
        park8_x(l,1) = park_position_struct(i).Lat;
        park8_y(l,1) = park_position_struct(i).Lon;
        l=l+1;
    elseif park_position_struct(i).Name == 9
        park9_x(m,1) = park_position_struct(i).Lat;
        park9_y(m,1) = park_position_struct(i).Lon;
        m=m+1;
    elseif park_position_struct(i).Name == 10
        park10_x(n,1) = park_position_struct(i).Lat;
        park10_y(n,1) = park_position_struct(i).Lon;
        n=n+1;
    elseif park_position_struct(i).Name == 11
        park11_x(o,1) = park_position_struct(i).Lat;
        park11_y(o,1) = park_position_struct(i).Lon;
        o=o+1;
    elseif park_position_struct(i).Name == 12
        park12_x(p,1) = park_position_struct(i).Lat;
        park12_y(p,1) = park_position_struct(i).Lon;
        p=p+1;
    elseif park_position_struct(i).Name == 13
        park13_x(q,1) = park_position_struct(i).Lat;
        park13_y(q,1) = park_position_struct(i).Lon;
        q=q+1;
    end
end

%Converting Lat & Lon coordinates into UTM
[park1_x_utm, park1_y_utm] = ll2utm([park1_x, park1_y]);
[park2_x_utm, park2_y_utm] = ll2utm([park2_x, park2_y]);
[park3_x_utm, park3_y_utm] = ll2utm([park3_x, park3_y]);
[park4_x_utm, park4_y_utm] = ll2utm([park4_x, park4_y]);
[park5_x_utm, park5_y_utm] = ll2utm([park5_x, park5_y]);
[park6_x_utm, park6_y_utm] = ll2utm([park6_x, park6_y]);
[park7_x_utm, park7_y_utm] = ll2utm([park7_x, park7_y]);
[park8_x_utm, park8_y_utm] = ll2utm([park8_x, park8_y]);
[park9_x_utm, park9_y_utm] = ll2utm([park9_x, park9_y]);
[park10_x_utm, park10_y_utm] = ll2utm([park10_x, park10_y]);
[park11_x_utm, park11_y_utm] = ll2utm([park11_x, park11_y]);
[park12_x_utm, park12_y_utm] = ll2utm([park12_x, park12_y]);
[park13_x_utm, park13_y_utm] = ll2utm([park13_x, park13_y]);

%shaping the sensitive places
park1 = polyshape([park1_x_utm],[park1_y_utm]);
park2 = polyshape([park2_x_utm],[park2_y_utm]);
park3 = polyshape([park3_x_utm],[park3_y_utm]);
park4 = polyshape([park4_x_utm],[park4_y_utm]);
park5 = polyshape([park5_x_utm],[park5_y_utm]);
park6 = polyshape([park6_x_utm],[park6_y_utm]);
park7 = polyshape([park7_x_utm],[park7_y_utm]);
park8 = polyshape([park8_x_utm],[park8_y_utm]);
park9 = polyshape([park9_x_utm],[park9_y_utm]);
park10 = polyshape([park10_x_utm],[park10_y_utm]);
park11 = polyshape([park11_x_utm],[park11_y_utm]);
park12 = polyshape([park12_x_utm],[park12_y_utm]);
park13 = polyshape([park13_x_utm],[park13_y_utm]);

parks = [park1; park2; park3; park4; park5; park6; park7; park8; park9; park10; park11; park12; park13;];

%build the sensitive area around each sensitive place (according to the
%Italian law, gNBs cannot be installed at a distance closer that 100m from
%each edge of each sensitive place
for i=1:length(parks)
    sensitive_areas(i,1) = polybuffer(parks(i), 100, 'JointType', 'miter');
end

clearvars -except parks torrino sensitive_areas d_min each_distance

run candidate_gNBs.m

%evaluate wheather a candidate Macro gNB site falls inside a sensistive area 
macro_in_sensitive_area = zeros(size(candidate_sites_macro,1), size(parks,1));
for i=1:size(parks,1)
    macro_in_sensitive_area(:,i) = inpolygon(candidate_sites_macro(:,1), candidate_sites_macro(:,2), sensitive_areas(i).Vertices(:,1), sensitive_areas(i).Vertices(:,2));
end

macro_in_sensitive_area = sum(macro_in_sensitive_area,2);
k=1;
to_delete = [];
for i=1:size(candidate_sites_macro,1)
    if macro_in_sensitive_area(i) ~= 0
        to_delete(k) = i;
        k=k+1;
    end
end

%delete from the list of Macro candidate sites the ones falling into at least one
%sensistive area
candidate_sites_macro(to_delete,:) = [];
nodes_macro = candidate_sites_macro;

%evaluate wheather a candidate micro gNB site falls inside a sensistive area 
micro_in_sensitive_area = zeros(size(candidate_sites_micro,1), size(parks,1));
for i=1:size(parks,1)
    micro_in_sensitive_area(:,i) = inpolygon(candidate_sites_micro(:,1), candidate_sites_micro(:,2), sensitive_areas(i).Vertices(:,1), sensitive_areas(i).Vertices(:,2));
end

micro_in_sensitive_area = sum(micro_in_sensitive_area,2);
k=1;
to_delete = [];
for i=1:size(candidate_sites_micro,1)
    if micro_in_sensitive_area(i) ~= 0
        to_delete(k) = i;
        k=k+1;
    end
end

%delete from the list of micro candidate sites the ones falling into at least one
%sensistive area
candidate_sites_micro(to_delete,:) = [];
nodes_micro = candidate_sites_micro;
