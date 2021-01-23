%Name of the file: micro_gNB.m
%Script version: v1.0
%Script authors: Luca Chiaraviglio, Cristian Di Paolo

%This script deals with the calculation of throughput and
%EMF parameters for each pixel under coverage of the considered micro
%gNBs

num_ofdm = 14; %number of ODFM symbols
tau = 3;  %OFDM symbols used for pilots
T_c = 500e-6;    %coherence time
cp_micro = 2.3e-6; %cyclic prefix
delta_F_micro = 30e3; %subcarrier spacing for micro cells
T_slot = (num_ofdm*(1/delta_F_micro)) + cp_micro; %slot duration
T_slot_micro=T_slot;
T_s = T_c/num_ofdm;  %symbol interval
T_s_micro=T_s;
T_pilot = tau*T_s;
T_pilot_micro=T_pilot;
T_u = 1/delta_F_micro; %useful symbol duration
T_u_micro=T_u;
T_g = T_s - T_u;  %guard interval
T_d = T_g; %Largest possible delay spread

%%parameters for power density computation per pixel related to micro cells
%Operating frequency: 3.7 GHz
Tx_gain_db = 15;
Tx_power_max_micro = 200;
alpha_stat = 0.25;
alpha_24 = 0.3;
Tx_power_effective = Tx_power_max_micro*alpha_stat*alpha_24;
Tx_power_effective_dbm = 10*log10(Tx_power_effective*1000);
Bw_h = 65; %horizontal BW
Bw_v = 10; %vertical BW
F = cosd(Bw_v)^2; %*sind(Bw_h)^2
Tx_gain = 10^(Tx_gain_db/10);
Tx_loss_db = 2.32;
EIRP_db = Tx_power_effective_dbm+Tx_gain_db-Tx_loss_db;
EIRP = 10.^(EIRP_db/10)/1000;
Z_0 = 377;

B_micro=80*10^6; %Total bandwidth
alpha_sec = 1; %number of sectors
alpha_sec_micro=alpha_sec;

%setting path loss variables for scenario 3GPP canyon UMi NLOS (TR 38.901 version 14)
gamma_micro = 3.19; %path loss exponent    
sigma_db_micro = 8.2; %std-dev variable for shadowing  

%setting the coverage distance for each micro gNB
cov_distance_micro = 200;

%log normal random variable generation and beta computation of the beta
%terms for the cell

rng(23595); %set the seed for the normal random number generator
distance = zeros(size(pixels_in_area,1), size(nodes_micro,1));
min_distance = zeros(size(pixels_in_area,1));
z_log = zeros(size(pixels_in_area, 1), size(extracted_site,2));
z = zeros(size(pixels_in_area, 1), size(extracted_site,2));
beta = zeros(size(pixels_in_area, 1), size(extracted_site,2));

%computing the distance between each pixel and the considered gNBs
for k=1:size(pixels_in_area,1)
    for j=1:size(nodes_micro,1)
        distance(k,j)=sqrt((sqrt((pixels_in_area(k,1)-nodes_micro(j,1))^2+(pixels_in_area(k,2)-nodes_micro(j,2))^2))^2+10^2);
        z_log(k,j) = normrnd(0,sigma_db_micro);
        z(k,j) = 10^(z_log(k,j)/10);
        beta(k,j)=z(k,j)/distance(k,j)^gamma_micro;
    end
    min_distance(k) = min(distance(k,:));
end

%assign each pixel to the closest gNB
coverage_matrix = zeros(size(pixels_in_area, 1), size(nodes_micro,1));
for i=1:size(distance,1)
    for j=1:size(nodes_micro,1)
        if distance(i,j) <= cov_distance_micro && min_distance(i) == distance(i,j) && sum(coverage_matrix(i,:),2) < 1
            coverage_matrix(i,j) = 1;
        else
            coverage_matrix(i,j) = 0;
        end
    end
end

%SIR & capacity computation
w=1;
SIR = [];
for k=1:size(pixels_in_area,1)
    beta_num = 0;
    sum_beta_den = 0;
    if sum(coverage_matrix(k,:)) == 1
        for j=1:size(nodes_micro,1)
            if (coverage_matrix(k,j) == 0)
                sum_beta_den = sum_beta_den+beta(k,j)^2;
            else
                beta_num = beta(k,j)^2;
            end
        end
        if(sum_beta_den > 0)
            SIR(k)=beta_num/sum_beta_den;
        else
            SIR(k)=beta_num/10^-32;
        end
    else
        SIR(k) = 0;
    end
end

w=1;
cov_SIR = [];
for i=1:length(SIR)
    if SIR(i) ~= 0
        cov_SIR(w) = SIR(i);
        w=w+1;
    end
end

SIR_db = 10*log10(SIR);
C_Mbps = [];
for k=1:size(SIR,2)
    C_Mbps(k) = ((B_micro/alpha_sec)*((T_slot-T_pilot)/T_slot)*(T_u/T_s)*log2(1+SIR(k)))/10^6;
end

w=1;
cov_C_Mbps = [];
for i=1:length(C_Mbps)
    if C_Mbps(i) ~= 0
        cov_C_Mbps(w) = C_Mbps(i);
        w=w+1;
    end
end

save('Throughput', 'SIR')

%EMF values computation
S_eq_mask = zeros(length(pixels_x), length(pixels_y), size(nodes_micro,1));
distance_pixel_mask = zeros(length(pixels_x), length(pixels_y), size(nodes_micro,1));

for i=1:length(pixels_x)
    for j=1:length(pixels_y)
        for k=1:size(nodes_micro,1)
            distance_pixel_mask(i,j,k)=sqrt((sqrt((pixels_x(i)-nodes_micro(k,1))^2+(pixels_y(j)-nodes_micro(k,2))^2))^2+10^2);
            S_eq_mask(i,j,k) = ((EIRP/(4*pi*distance_pixel_mask(i,j,k)^2))*F);
        end
    end
end

S_eq_micro = zeros(size(pixel_coordinate,1), size(nodes_micro,1));
distance_pixel = zeros(size(pixel_coordinate,1), size(nodes_micro,1));

for i=1:size(pixel_coordinate,1)
    for k=1:size(nodes_micro,1)
        distance_pixel(i,k)=sqrt((sqrt((pixel_coordinate(i,1)-nodes_micro(k,1))^2+(pixel_coordinate(i,2)-nodes_micro(k,2))^2))^2+10^2);
        S_eq_micro(i,k) = ((EIRP/(4*pi*distance_pixel(i,k)^2))*F);
    end
end

S_eq_mask_micro = sum(S_eq_mask, 3);
E_mask_micro = sqrt(S_eq_mask_micro.*Z_0);

exclusion_zone = zeros(size(E_mask_micro,1), size(E_mask_micro,2));
for i=1:size(E_mask_micro,1)
    for j=1:size(E_mask_micro,2)
        if (E_mask_micro(i,j) <= 6)
            exclusion_zone(i,j) = 0;
        else
            exclusion_zone(i,j) = 1;
        end
    end
end

%saving useful variables to be used later 
save('EMF variables', 'S_eq_mask_micro', 'S_eq_micro')