%Name of the file: PLATEA_main.m
%Script version: v1.0
%Script authors: Luca Chiaraviglio, Cristian Di Paolo

%%This script is designed to evaluate the best scenario
%%planning for the Torrino-Mezzocammino neighborhood out of a range of
%%possible combinations. This is done by exploiting all the potentials
%%of the objectvie function in terms of minimizing the total
%%scenario costs and maximizing the coverage, always complying with the 
%%EMF restrictions and regulations 

close all
clear all
clc

tic;

%generate the map of the considered area and divide it in 10x10m pixels
run locate_sensitive_places.m

pixels_x = [2.8525*10^5:10:2.8775*10^5];
pixels_y = [4.62925*10^6:10:4.63175*10^6];

load('pixel_coordinates.mat');

inside_torrino = inpolygon(pixel_coordinate(:,1), pixel_coordinate(:,2), torrino.Vertices(:,1), torrino.Vertices(:,2));

k=1;
w=1;
for i=1:length(inside_torrino)
    if inside_torrino(i) == 0
        pixels_out(k,:) = pixel_coordinate(i,:);
        k=k+1;
    else
        pixels_in_area(w,:) = pixel_coordinate(i,:);
        w=w+1;
    end
end

%Starting the iteration over the micro cells
curr_scenario = 1; %varialble for controlling obj_min
for n_micro=1:size(candidate_sites_micro,1)-1
    
    %counting micro to save the execution times
    n_micro_vector(n_micro) = n_micro;

    %setting scenarios fo executing permutations
    extracted_site = [];
    extracted_site(1,:) = randperm(size(candidate_sites_micro,1),n_micro);
    
    for n_permutations=2:size(candidate_sites_micro,1)-1
        extracted = randperm(size(candidate_sites_micro,1),n_micro);
        if ~ismember(extracted_site,extracted,'rows')
            extracted_site(n_permutations,:) = extracted;
        else
            extracted_site(n_permutations,:) = randperm(size(candidate_sites_micro,1),n_micro);
        end
    end
    
    nodes_micro = zeros(size(extracted_site,2),2);

    %executing permutations
    for n_scenario=1:size(extracted_site,1)
        for b=1:size(extracted_site,2)
            nodes_micro(b,:) = [candidate_sites_micro(extracted_site(n_scenario,b),:);]; 
        end
        
        %nodes position matrix
        area_exclusion = zeros(size(pixels_x,2), size(pixels_y,2));
        distance_micro = [];
        for k=1:size(nodes_micro,1)
            for i=1:length(pixels_x)
                for j=1:length(pixels_y)
                    distance_micro(i,j,k) = sqrt((pixels_x(i)-nodes_micro(k,1))^2+(pixels_y(j)-nodes_micro(k,2))^2);
                end
            end
        end
        
        %setting the exclusion zones
        for i=1:size(nodes_micro,1)
            for j=1:length(pixels_x)
                for k=1:length(pixels_y)
                    if distance_micro(j,k,i) <= 11.5  
                        area_exclusion(j,k) = 1;
                    end
                end
            end
        end
       
        run micro_gNB.m

        %setting control vartiables
        y_lf = ones(1,size(nodes_micro,1));
        x_plf = zeros(size(pixels_in_area,1), size(nodes_micro,1));
        m_l = zeros(1,size(nodes_micro,1));
        N_ser = 2;
        x_plf(:,:,1) = coverage_matrix;

        %%check that the distance between the pixel and the BS is lower than the
        %%maximum coverage of the BS itself
        beyond_coverage = [];
        k=1;
        for i=1:size(distance,1)
            if sum(x_plf(i,:,1)) == 0
                beyond_coverage(k,:) = pixel_coordinate(i,:);
                k=k+1;
            end
        end

        %%check that each pixel can be served by up to N_ser BSs
        for i=1:size(x_plf,1)
            if sum(x_plf(i,:,1))>N_ser
                continue
            end
        end

        %%check if the SIR value for each pixel is enough to guarantee a min Capacity of
        %%30Mbit/s
        C_min = 30*10^6;
        SIR_min = 2^(C_min/((B_micro/alpha_sec)*((T_slot-T_pilot)/T_slot)*(T_u/T_s)))-1;
        below_SIR_before_macro = [];

        k = 1;

        for i=1:size(x_plf,1)
            if sum(x_plf(i,:,1)) == 1
                if SIR(i) < SIR_min
                    below_SIR_before_macro(k,:) = pixel_coordinate(i,:);
                    x_plf(i,:,1) = 0;
                    k=k+1;
                end
            end
        end

        %%check if a pixel falls in the exclusion zone of an installed BS
        k=1;
        w_p = [];
        for i=1:size(area_exclusion,1)
            for j=1:size(area_exclusion,2)
                w_p(k) = area_exclusion(i,j);
                k = k+1;
            end
        end

        k=1;
        exclusion_vector = [];
        for i=1:size(exclusion_zone,1)
            for j=1:size(exclusion_zone,2)
                exclusion_vector(k) = exclusion_zone(i,j);
                k = k+1;
            end
        end
        
        k=1;
        inside_exclusion_zone = [];
        for i=1:length(exclusion_vector)
            if exclusion_vector(i) == 1
                inside_exclusion_zone(k,:) = pixel_coordinate(i,:);
                x_plf(i,:,1) = 0;
                k=k+1;
            end
        end

        %%check on the power density limits
        L_res = 0.1; %power density limit over residential areas
        outside_exclusion = true;
        for i=1:length(w_p)
            if ((1-w_p(i))*sum(sum(S_eq_micro(i,:,:),2)))/L_res > 1
                outside_exclusion = false;
            end
        end   
        
        %perform the objective functions calculation only if the the area
        %under coverage does not include pixels inside the exclusion zone
        if outside_exclusion
            
            %%computation of total cost for each installed site
            installation_costs = [14852, 20101]; %installation costs for micro and macro, respectively
            equipment_cost = [2791, 45673];  %equipment costs for micro and macro, respectively
            
            C_lf_SITE = [];
            C_f_EQUIP = [];
            for i=1:size(y_lf,1)
                C_lf_SITE(i,:) = installation_costs(i).*y_lf(i,:); %site installation cost at location l, depending on the kond of cell
                C_f_EQUIP(i,:) = equipment_cost(i).*y_lf(i,:); %equipment cost of a BS installed depending on type
            end
            
            sum_C_lf_SITE = sum(C_lf_SITE,1);
            sum_C_f_EQUIP = sum(C_f_EQUIP,1);
            
            for i=1:size(nodes_micro,1)
                m_l(i) = sum_C_lf_SITE(i) + sum_C_f_EQUIP(i); 
            end
            
            %%objective function calculation
            alpha_micro = 2; %weight factor
            obj(curr_scenario,n_scenario) = sum(m_l,2)-alpha_micro*sum(sum(sum(x_plf,3),1),2); %%setting alpha to get obj_min
        else 
            obj(curr_scenario,n_scenario) = NaN;
        end
    end
    if isnan(min(obj(curr_scenario,:))) 
        continue
    else
        [obj_min(curr_scenario), index_min(curr_scenario)] = min(obj(curr_scenario,:));
        
        %saving the best scenario
        best_scenario(curr_scenario,:) = {candidate_sites_micro(extracted_site(index_min(curr_scenario),:),:)};
        curr_scenario = curr_scenario+1;
    end
end

best_scenario(curr_scenario,:) = {candidate_sites_micro};

%reset obj variables for overall analysis
obj = [];
obj_min = [];
index_min = [];

%defining variables for summarizing the best scenario relevant
%characteristics
var_1 = 'Number_of_micro_cells_installed';
var_2 = 'Number_of_macro_cells_installed';
var_3 = 'Number_of_pixels_served_by_macro_cells';
var_4 = 'Number_of_pixels_served_by_micro_cells';
var_5 = 'alpha_value_for_micro';
var_6 = 'alpha_value_for_macro';
var_7 = 'Objective_function_value';
var_8 = 'Total_scenario_cost';
var_9 = 'BSs_coordinates';
var_10 = 'Average_throughput_per_covered_user_Mbps';
var_11 = 'Micro_costs';
var_12 = 'Macro_costs';
var_13 = 'Percentage_pixels_not_served';
var_14 = 'S_eq_avg_over_L_res';
var_15 = 'Average_throughput_Mbps';
var_16 = 'Average_EMF_value';
var_17 = 'Average_througput_micro_cov';
var_18 = 'Average_througput_macro_cov';

%set a vector of costs for a single micro or macro
costs = installation_costs + equipment_cost;

%adding Macro coverage to each best sceanario just obtained
for each_best=1:size(best_scenario,1)
    
    n_macro_vector(each_best) = each_best;

    %assigning best micro scenario for each iteration
    nodes_micro = best_scenario{each_best,1};
    
    run micro_gNB.m
    
    %setting control variables
    y_lf = ones(1,size(nodes_micro,1));
    x_plf = zeros(size(pixels_in_area,1), size(nodes_micro,1));
    m_l = zeros(1,size(nodes_micro,1));
    N_ser = 2;
    x_plf(:,:,1) = coverage_matrix;
    
    %%check that the distance between the pixel and the BS is lower than the
    %%maximum coverage of the BS itself
    beyond_coverage = [];
    k=1;
    for i=1:size(distance,1)
        if sum(x_plf(i,:,1)) == 0
            beyond_coverage(k,:) = pixel_coordinate(i,:);
            k=k+1;
        end
    end
    
    below_SIR_before_macro = [];
    k=1;
    for i=1:size(x_plf,1)
        if sum(x_plf(i,:,1)) == 1
            if SIR(i) < SIR_min
                below_SIR_before_macro(k,:) = pixel_coordinate(i,:);
                x_plf(i,:,1) = 0;
                k=k+1;
            end
        end
    end
    
    save('x_plf', 'x_plf')
    
    k=1;
    flag = true;
    n_macro = 0;
    while flag
        n_macro = n_macro + 1;
        
        extracted_site = [];
        extracted_site(1,:) = randperm(size(candidate_sites_macro,1),n_macro);
        
        if n_macro < size(candidate_sites_macro,1)
            for n_permutations=2:size(candidate_sites_macro,1)
                extracted = randperm(size(candidate_sites_macro,1),n_macro);
                if ~ismember(extracted_site,extracted,'rows')
                    extracted_site(n_permutations,:) = extracted;
                else
                    extracted_site(n_permutations,:) = randperm(size(candidate_sites_macro,1),n_macro);
                end
            end
        else
            extracted_site = randperm(size(candidate_sites_macro,1),n_macro);
        end
        
        nodes_macro = zeros(size(extracted_site,2),2);
        
        for n_scenario=1:size(extracted_site,1)
            for b=1:size(extracted_site,2)
                nodes_macro(b,:) = [candidate_sites_macro(extracted_site(n_scenario,b),:);]; %candidate_sites(extracted_site(a,2),:); candidate_sites(extracted_site(a,3),:); candidate_sites(extracted_site(a,4),:);];
            end
            run macro_gNB.m
            load('New_served.mat')
            
            %nodes position matrix
            area_exclusion = zeros(size(pixels_x,2), size(pixels_y,2));
            distance_micro = [];
            for k=1:size(nodes_micro,1)
                for i=1:length(pixels_x)
                    for j=1:length(pixels_y)
                        distance_micro(i,j,k) = sqrt((pixels_x(i)-nodes_micro(k,1))^2+(pixels_y(j)-nodes_micro(k,2))^2);
                    end
                end
            end
            
            for i=1:size(nodes_micro,1)
                for j=1:length(pixels_x)
                    for k=1:length(pixels_y)
                        if distance_micro(j,k,i) <= 12
                            area_exclusion(j,k) = 1;
                        end
                    end
                end
            end
            
            distance_macro = [];
            for k=1:size(nodes_macro,1)
                for i=1:length(pixels_x)
                    for j=1:length(pixels_y)
                        distance_macro(i,j,k) = sqrt((pixels_x(i)-nodes_macro(k,1))^2+(pixels_y(j)-nodes_macro(k,2))^2);
                    end
                end
            end
            
            for i=1:size(nodes_macro,1)
                for j=1:length(pixels_x)
                    for k=1:length(pixels_y)
                        if distance_macro(j,k,i) <= 5
                            area_exclusion(j,k) = 1;
                        end
                    end
                end
            end
            
            %saving the average throughput value per each scenario
            mean_throughput_covered(n_macro,n_scenario,each_best) = mean(cov_C_Mbps);
            mean_throughput_overall(n_macro,n_scenario,each_best) = mean(C_Mbps);
            mean_throughput_micro(n_macro,n_scenario,each_best) = mean(cov_micro_C_Mbps);
            mean_throughput_macro(n_macro,n_scenario,each_best) = mean(cov_macro_C_Mbps);
            
            %update pixels not served, taking into account
            %both kinds of BSs (at the end not_served already includes
            %also pixels covered by micro but under SIR min)
            not_served = [];
            k = 1;
            for i=1:size(x_plf,1)
                if sum(x_plf(i,:,1)) == 0 && sum(x_plf(i,:,2)) == 0
                    not_served(k,:) = pixels_in_area(i,:);
                    k=k+1;
                end
            end
            
            if size(not_served,1) == 0
                flag = false;
            end
            
            %counting how many pixels are not served for each iteration
            n_not_served(n_macro,n_scenario,each_best) = size(not_served,1);
            
            %differentiating pixels served by micro and Macro
            k=1;
            w=1;
            served_by_micro = [];
            served_by_Macro = [];
            for i=1:size(x_plf,1)
                if sum(x_plf(i,:,1)) ~= 0
                    served_by_micro(k,:) = pixels_in_area(i,:);
                    k=k+1;
                elseif sum(x_plf(i,:,2)) ~= 0
                    served_by_Macro(w,:) = pixels_in_area(i,:);
                    w=w+1;
                end
            end
            
            %counting how many pixels are served by micro and Macro BSs
            n_served_by_micro(n_macro,n_scenario,each_best) = size(served_by_micro,1);
            n_served_by_Macro(n_macro,n_scenario,each_best) = size(served_by_Macro,1);
            
            %updating the list of pixels that remain below SIR min also
            %after the installation of Macro
            below_SIR_after_macro = below_SIR_before_macro;
            if size(below_SIR_before_macro,1) > 0
                below = [];
                below = ismember(below_SIR_before_macro,served_by_Macro, 'rows');
                k=1;
                to_delete = [];
                if sum(below,1) > 0
                    for i=1:size(below,1)
                        if below(i) == 1
                            to_delete(k) = i;
                            k=k+1;
                        end
                    end
                    below_SIR_after_macro(to_delete,:) = [];
                end
            end
            
            %%check that each pixel can be served by up to N_ser BSs
            for i=1:size(x_plf,1)
                if sum(x_plf(i,:,1)) + sum(x_plf(i,:,2)) > N_ser
                    warning('Pixel %d is served by more than 2 BSs', i);
                end
            end
            
            y_lf = [zeros(1,size(nodes_macro,1)) ones(1,size(nodes_micro,1)); ones(1,size(nodes_macro,1)) zeros(1,size(nodes_micro,1))];
            u_l = ones(1,size(nodes,1));
            
            %EMF values assessment
            S_tot_mask = zeros(size(E_mask_micro,1), size(E_mask_micro,2));
            
            load('EMF variables')
            for a=1:size(S_eq_mask_macro)
                for b=1:size(S_eq_mask_macro)
                    S_tot_mask(a,b) = S_eq_mask_micro(a,b) + S_eq_mask_macro(a,b);
                end
            end
            
            E_tot_mask = sqrt(S_tot_mask.*Z_0);
            
            %computing average EMF value over the considered scenario
            E_avg(n_macro,n_scenario,each_best) = mean(E_tot_mask,'all');
            
            L_res = 0.1; %power density limit over residential areas
            S_avg(n_macro,n_scenario,each_best) = mean(S_tot_mask/L_res, 'all');
            
            S_eq_plf = zeros(size(S_eq_macro,1), size(nodes,1), size(y_lf,1));
            S_eq_plf(:,:,1) = [zeros(size(S_eq_macro,1),size(S_eq_macro,2)) S_eq_micro];
            S_eq_plf(:,:,2) = [S_eq_macro zeros(size(S_eq_micro,1),size(S_eq_micro,2))];
            
            exclusion_zone = zeros(size(E_tot_mask,1),size(E_tot_mask,1));
            for i=1:size(E_tot_mask,1)
                for j=1:size(E_tot_mask,2)
                    if (E_tot_mask(i,j) <= 6)
                        exclusion_zone(i,j) = 0;
                    else
                        exclusion_zone(i,j) = 1;
                    end
                end
            end
            
            %%check if a pixel falls in the exclusion zone of an installed BS
            k=1;
            w_p = [];
            for i=1:size(area_exclusion,1)
                for j=1:size(area_exclusion,2)
                    w_p(k) = area_exclusion(i,j);
                    k = k+1;
                end
            end
            
            k=1;
            exclusion_vector = [];
            for i=1:size(exclusion_zone,1)
                for j=1:size(exclusion_zone,2)
                    exclusion_vector(k) = exclusion_zone(i,j);
                    k = k+1;
                end
            end
            
            k=1;
            inside_exclusion_zone = [];
            for i=1:length(exclusion_vector)
                if exclusion_vector(i) == 1
                    inside_exclusion_zone(k,:) = pixel_coordinate(i,:);
                    k=k+1;
                end
            end
            
            %%check on the power density limits
            outside_exclusion = true;
            for i=1:length(w_p)
                if ((1-w_p(i))*sum(sum(S_eq_plf(i,:,:),2)))/L_res > 1
                    outside_exclusion = false;
                end
            end
            
            C_lf_SITE = [];
            C_f_EQUIP = [];
            for i=1:size(y_lf,1)
                C_lf_SITE(i,:) = installation_costs(i).*y_lf(i,:); %site installation cost at location l, depending on the kond of cell
                C_f_EQUIP(i,:) = equipment_cost(i).*y_lf(i,:); %equipment cost of a BS installed depending on type
            end
            
            sum_C_lf_SITE = sum(C_lf_SITE,1);
            sum_C_f_EQUIP = sum(C_f_EQUIP,1);
            
            for i=1:size(nodes,1)
                m_l(i) = sum_C_lf_SITE(i) + sum_C_f_EQUIP(i); % + tower_cost
            end
            
            total_cost(n_macro, each_best) = sum(m_l,2);
            
            %%overall objective function calculation
            alpha_micro = logspace(1,7,200); %weight factor for micro
            alpha_macro = logspace(1,4,200); %weight factor for mccro
            k=1;
            for n_alpha_micro=1:length(alpha_micro)
                for n_alpha_macro=1:length(alpha_macro)
                    obj(n_macro,n_scenario,each_best,k) = sum(m_l,2)-(alpha_micro(n_alpha_micro)*sum(sum(x_plf(:,:,1))) + alpha_macro(n_alpha_macro)*sum(sum(x_plf(:,:,2)))); %%setting alpha to get obj_min
                    k=k+1;
                end
            end
        end
        if n_macro < size(candidate_sites_macro,1)
            k=1;
            for n_alpha_micro=1:length(alpha_micro)
                for n_alpha_macro=1:length(alpha_macro)
                    [obj_min(each_best, n_macro, k), index_min(each_best, n_macro, k)] = min(obj(n_macro,:,each_best,k));
                    
                    %saving the most relevant variables for each overall best
                    %scenario
                    summary(each_best, n_macro, k) = struct(var_1, size(best_scenario{each_best},1), var_2, n_macro, var_3, n_served_by_Macro(n_macro, index_min(each_best, n_macro, k), each_best), var_4, n_served_by_micro(n_macro, index_min(each_best, n_macro, k), each_best), var_5, alpha_micro(n_alpha_micro), var_6, alpha_macro(n_alpha_macro), var_7, obj_min(each_best, n_macro, k), var_8, total_cost(n_macro, each_best), var_9, {vertcat(best_scenario{each_best,1}, candidate_sites_macro(extracted_site(index_min(each_best, n_macro),:),:))}, var_10, mean_throughput_covered(n_macro, index_min(each_best, n_macro, k), each_best), var_11, costs(1,1)*size(best_scenario{each_best},1), var_12, costs(1,2)*n_macro, var_13, n_not_served(n_macro, index_min(each_best, n_macro, k), each_best)/size(pixels_in_area,1)*100, var_14, S_avg(n_macro, index_min(each_best, n_macro, k), each_best), var_15, mean_throughput_overall(n_macro, index_min(each_best, n_macro, k), each_best), var_16, E_avg(n_macro, index_min(each_best, n_macro, k), each_best), var_17, mean_throughput_micro(n_macro, index_min(each_best, n_macro), each_best), var_18, mean_throughput_macro(n_macro, index_min(each_best, n_macro), each_best));
                    k=k+1;
                end
            end
        else
            k=1;
            for n_alpha_micro=1:length(alpha_micro)
                for n_alpha_macro=1:length(alpha_macro)
                    obj_min(each_best, n_macro, k) = obj(n_macro,n_scenario,each_best,k);
                    index_min(each_best, n_macro, k) = 1;
                    summary(each_best, n_macro, k) = struct(var_1, size(best_scenario{each_best},1), var_2, n_macro, var_3, n_served_by_Macro(n_macro, index_min(each_best, n_macro, k), each_best), var_4, n_served_by_micro(n_macro, index_min(each_best, n_macro, k), each_best), var_5, alpha_micro(n_alpha_micro), var_6, alpha_macro(n_alpha_macro), var_7, obj_min(each_best, n_macro, k), var_8, total_cost(n_macro, each_best), var_9, {vertcat(best_scenario{each_best,1}, candidate_sites_macro(extracted_site(index_min(each_best, n_macro),:),:))}, var_10, mean_throughput_covered(n_macro, index_min(each_best, n_macro, k), each_best), var_11, costs(1,1)*size(best_scenario{each_best},1), var_12, costs(1,2)*n_macro, var_13, n_not_served(n_macro, index_min(each_best, n_macro, k), each_best)/size(pixels_in_area,1)*100, var_14, S_avg(n_macro, index_min(each_best, n_macro, k), each_best), var_15, mean_throughput_overall(n_macro, index_min(each_best, n_macro, k), each_best), var_16, E_avg(n_macro, index_min(each_best, n_macro, k), each_best), var_17, mean_throughput_micro(n_macro, index_min(each_best, n_macro), each_best), var_18, mean_throughput_macro(n_macro, index_min(each_best, n_macro), each_best));
                    k=k+1;
                end
            end
        end
    end
end

%%saving only overall scenarios with the min value of alpha
for i=1:size(summary,1)
    for j=1:size(summary,2)
        for k=1:size(summary,3)
            if summary(i,j,k).Objective_function_value == min(obj_min(:,:,k),[],'all')
                summary_alpha(k) = summary(i,j,k);
            end
        end
    end
end

%save the struct containing the elements needed to be further analyzed
%running the alpha_heatmaps.m file
save('Summary','summary_alpha', 'alpha_macro', 'alpha_micro')

toc;
