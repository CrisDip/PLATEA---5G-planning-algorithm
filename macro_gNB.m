%Name of the file: macro_gNB.m
%Script version: v1.0
%Script authors: Luca Chiaraviglio, Cristian Di Paolo

%This script deals with the calculation of throughput and
%EMF parameters for each pixel under coverage of the considered Macro
%gNBs

num_ofdm = 14; %number of ODFM symbols
tau = 3;  %OFDM symbols used for pilots
T_c = 500e-6;    %coherence time
cp_macro = 4.7e-6; %cyclic prefix
delta_F_macro = 15e3; %subcarrier spacing for micro cells
T_slot = (num_ofdm*(1/delta_F_macro)) + cp_macro; %slot duration
T_s = T_c/num_ofdm;  %symbol interval
T_pilot = tau*T_s;
T_u = 1/delta_F_macro; %useful symbol duration
T_g = T_s - T_u;  %guard interval
T_d = T_g; %Largest possible delay spread

%Operating frequency: 700 MHz
B_macro=20*10^6; %Total bandwidth
alpha_sec=1; %frequency reuse factor

%setting path loss variables for scenario 3GPP canyon UMa NLOS (TR 38.901 version 14)
gamma_macro = 3; %path loss exponent    
sigma_db_macro = 6.8; %std-dev variable for shadowing  

%setting the coordinates for each pixel
load('x_plf')

%checking which pixels are not served
not_served = [];
for i=1:size(x_plf,1)
    if sum(x_plf(i,:,1)) == 0
        not_served(i,:) = pixels_in_area(i,:);
    end
end

%setting the coverage distance for each Macro gNB
cov_distance_macro = 900;

%log normal random variable generation and beta computation of the beta
%terms for the cell

rng(101197); %set the seed for the normal random number generator
distance = zeros(size(pixels_in_area,1), size(nodes_macro,1));
z_log = zeros(size(pixels_in_area, 1), size(extracted_site,2));
z = zeros(size(pixels_in_area, 1), size(extracted_site,2));
beta = zeros(size(pixels_in_area, 1), size(extracted_site,2));
min_distance = zeros(size(pixels_in_area,1));

%computing the distance between each pixel and the considered gNBs
for k=1:size(pixels_in_area,1)
    for j=1:size(nodes_macro,1)
        distance(k,j)=sqrt((sqrt((pixels_in_area(k,1)-nodes_macro(j,1))^2+(pixels_in_area(k,2)-nodes_macro(j,2))^2))^2+25^2);
        z_log(k,j) = normrnd(0,sigma_db_macro);
        z(k,j) = 10^(z_log(k,j)/10);
        beta(k,j)=z(k,j)/distance(k,j)^gamma_macro;
    end
    min_distance(k) = min(distance(k,:));
end

%assign each pixel to the closest BS
coverage_matrix = zeros(size(pixels_in_area, 1), size(nodes_macro,1));
for i=1:size(distance,1)
    for j=1:size(nodes_macro,1)
        if distance(i,j) <= cov_distance_macro && distance(i,j) ~= 0 && min_distance(i) == distance(i,j) && sum(coverage_matrix(i,:),2) < 1
            coverage_matrix(i,j) = 1;
        else
            coverage_matrix(i,j) = 0;
        end
    end
end

%updating the nodes and pixels considered, after Macro gNBs installation
nodes=[nodes_macro; nodes_micro];

x_plf = [x_plf zeros(size(x_plf,1),size(nodes_macro,1))];

k=1;
new_served = [];
for j=1:size(nodes_macro,1)
    for i=1:length(coverage_matrix)
        if coverage_matrix(i,j) == 1
            new_served(k,:) = pixels_in_area(i,:);
            x_plf(i,size(nodes_micro,1)+j,2) = 1;
            k=k+1;
        end
    end
end
save('New_served', 'new_served', 'x_plf')

%SIR computation
load('Throughput')
for k=1:size(pixels_in_area,1)
    beta_num = 0;
    sum_beta_den = 0;
    if sum(x_plf(k,:,2)) == 1 
        for j=1:size(nodes_macro,1)
            if (x_plf(k,size(nodes_micro,1)+j,2) == 0)
                sum_beta_den = sum_beta_den+beta(k,j)^2;
            else
                beta_num = beta(k,j)^2;
            end
        end
        if(sum_beta_den > 0)
            SIR(2,k)=beta_num/sum_beta_den;
        else
            SIR(2,k)=beta_num/10^-32;
        end
    elseif sum(x_plf(k,:,2)) == 0 && sum(x_plf(k,:,1)) == 0
        SIR(2,k) = 0;
    end
end

SIR_tot = sum(SIR,1);
w=1;
cov_SIR = [];
for i=1:length(SIR_tot)
    if SIR_tot(i) ~= 0
        cov_SIR(w) = SIR_tot(i);
        w=w+1;
    end
end

SIR_db = [];
SIR_db = 10*log10(cov_SIR);

C_Mbps = [];
for k=1:size(SIR,2)
    C_Mbps(1,k) = ((B_micro/alpha_sec_micro)*((T_slot_micro-T_pilot_micro)/T_slot_micro)*(T_u_micro/T_s_micro)*log2(1+SIR(1,k)))/10^6;
    C_Mbps(2,k) = ((B_macro/alpha_sec)*((T_slot-T_pilot)/T_slot)*(T_u/T_s)*log2(1+SIR(2,k)))/10^6;
end

w=1;
cov_micro_C_Mbps = [];
for j=1:size(C_Mbps,2)
    if C_Mbps(1,j) ~= 0
        cov_micro_C_Mbps(w) = C_Mbps(1,j);
        w=w+1;
    end
end

w=1;
cov_macro_C_Mbps = [];
for j=1:size(C_Mbps,2)
    if C_Mbps(2,j) ~= 0
        cov_macro_C_Mbps(w) = C_Mbps(2,j);
        w=w+1;
    end
end

C_Mbps = sum(C_Mbps,1);

w=1;
cov_C_Mbps = [];
for i=1:length(C_Mbps)
    if C_Mbps(i) ~= 0
        cov_C_Mbps(w) = C_Mbps(i);
        w=w+1;
    end
end

%%power density computation per pixel
Tx_gain_db = 15;
Tx_power_max_macro = 65;
alpha_24 = 0.3;
Tx_power_effective = Tx_power_max_macro*alpha_24;
Tx_power_effective_dbm = 10*log10(Tx_power_effective*1000);
Bw_h = 65; %horizontal BW
Bw_v = 10; %vertical BW
F = cosd(Bw_v)^2; %*sind(Bw_h)^2
Tx_gain = 10^(Tx_gain_db/10);
Tx_loss_db = 2.32;
EIRP_db = Tx_power_effective_dbm+Tx_gain_db-Tx_loss_db;
EIRP = 10.^(EIRP_db/10)/1000;
Z_0 = 377;

S_eq_mask = zeros(length(pixels_x), length(pixels_y), size(nodes_macro,1));
distance_pixel_mask = zeros(length(pixels_x), length(pixels_y), size(nodes_macro,1));

for i=1:length(pixels_x)
    for j=1:length(pixels_y)
        for k=1:size(nodes_macro,1)
            distance_pixel_mask(i,j,k)=sqrt((sqrt((pixels_x(i)-nodes_macro(k,1))^2+(pixels_y(j)-nodes_macro(k,2))^2))^2+25^2);
            S_eq_mask(i,j,k) = ((EIRP/(4*pi*distance_pixel_mask(i,j,k)^2))*F);
        end
    end
end

S_eq_macro = zeros(size(pixel_coordinate,1), size(nodes_macro,1));
distance_pixel = zeros(size(pixel_coordinate,1), size(nodes_macro,1));

for i=1:size(pixel_coordinate,1)
    for k=1:size(nodes_macro,1)
        distance_pixel(i,k)=sqrt((sqrt((pixel_coordinate(i,1)-nodes_macro(k,1))^2+(pixel_coordinate(i,2)-nodes_macro(k,2))^2))^2+25^2);
        S_eq_macro(i,k) = ((EIRP/(4*pi*distance_pixel(i,k)^2))*F);
    end
end

S_eq_mask_macro = sum(S_eq_mask, 3);