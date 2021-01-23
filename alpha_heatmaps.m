%Name of the file: alpha_heatmaps.m
%Script version: v1.0
%Script authors: Luca Chiaraviglio, Cristian Di Paolo

%%This script is designed to further analize the
%%relationship between the main variables involved in the planning, such as
%%capacity, total costs and average EMF values, and the variation of the 
%%alpha weight parameter.

close all
clear all

%load the struct in which all the information has been saved
load('Summary')

%show the relationship between micro insatllation costs variation and alpha variation
micro_costs_matrix = reshape([summary_alpha.Micro_costs], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1 = axes;
imagesc(alpha_macro, alpha_micro, micro_costs_matrix)
set(h1,'YDir','normal')
set(h1,'YScale', 'log')
set(h1,'XScale', 'log')
title('Micro cells costs vs alpha')
h = colorbar;
ylabel(h, 'Euros')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between Macro installation costs variation and alpha variation
macro_costs_matrix = reshape([summary_alpha.Macro_costs], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1 = axes;
imagesc(alpha_macro, alpha_micro, (macro_costs_matrix.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Macro cells costs vs alpha')
h = colorbar;
ylabel(h, 'Euros')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between obtained throughput variation and alpha
%variation (only covered pixels)
throughput_matrix_covered = reshape([summary_alpha.Average_throughput_per_covered_user_Mbps], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (throughput_matrix_covered.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Average throughput (only covered) vs alpha')
h = colorbar;
ylabel(h, 'Mbps')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between obtained throughput variation and alpha
%variation (all pixels)
throughput_matrix_overall = reshape([summary_alpha.Average_throughput_Mbps], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (throughput_matrix_overall.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Average throughput (all pixels) vs alpha')
h = colorbar;
ylabel(h, 'Mbps')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between overall installation costs variation and alpha variation
total_costs_matrix = reshape([summary_alpha.Total_scenario_cost], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (total_costs_matrix.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Macro + micro cells costs vs alpha')
h = colorbar;
ylabel(h, 'Euros')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between the number of micro installed and alpha variation
n_micro_matrix = reshape([summary_alpha.Number_of_micro_cells_installed], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, n_micro_matrix)
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Number of micro installed vs alpha')
h = colorbar;
ylabel(h, 'N.')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between the number of Macro installed and alpha variation
n_macro_matrix = reshape([summary_alpha.Number_of_macro_cells_installed], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (n_macro_matrix.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Number of Macro installed vs alpha')
h = colorbar;
ylabel(h, 'N.')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between the number of number of pixels not served and alpha variation
not_served_matrix = reshape([summary_alpha.Percentage_pixels_not_served], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (not_served_matrix.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Percentage of pixels not served vs alpha')
h = colorbar;
ylabel(h, '%')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between avg power density and alpha variation
S_eq_avg_matrix = reshape([summary_alpha.S_eq_avg_over_L_res], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (S_eq_avg_matrix.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('Power density average vs alpha')
h = colorbar;
ylabel(h, '%')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%show the relationship between the avg EMF values and alpha variation
EMF_average = reshape([summary_alpha.Average_EMF_value], size(alpha_macro,2), size(alpha_micro,2));

figure()
h1=axes;
imagesc(alpha_macro, alpha_micro, (EMF_average.').')
set(h1,'YDir','normal', 'YScale', 'log')
set(h1,'XScale', 'log')
title('EMF average')
h = colorbar;
ylabel(h, '[V/m]')
xlabel('\alpha_{24}');
ylabel('\alpha_{stat}');

%plot the throughput evolution according to the number of micro installed
k=1;
for i=1:size(summary_alpha,2)
    if summary_alpha(i).Number_of_macro_cells_installed == 3 && summary_alpha(i).Number_of_micro_cells_installed ~= 1
        micro(k) = summary_alpha(i).Number_of_micro_cells_installed;
        capacity(k) = summary_alpha(i).Average_throughput_per_covered_user_Mbps;
        k=k+1;
    end
end

% micro = unique(micro);
% capacity = unique(capacity);

figure()
h_3 = plot(micro, capacity);
set(h_3, 'LineWidth',1.5)
grid on
xlabel('n. of micro installed')
ylabel('Capacity [Mbps]')
title('Evaluation of avg capacity. Totally covered scenario')
    