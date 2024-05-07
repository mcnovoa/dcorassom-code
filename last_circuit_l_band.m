clc;clear all;close all;

T = 20:0.1:32;
SSS = 0:0.1:45;
j = 1;


vout_r = [];
pin_r = [];
result = [0,0,0];
tb_h = [];
ta = [];
TK = [];
ta_n = [];
p_tbh = [];

vout_a = [];
pin_a = [];

for i = 1:length(SSS)
    for j = 1:length(T)
        result = volt_function(28, SSS(i));
        vout_r(i) = result(2);
        pin_r(i) = result(1);
        tb_h(i) = result(3);
        ta(i) = result(4);
        vout_a(i) = result(5);
        pin_a(i) = result(6);
        ta_n(i) = result(7);
        p_tbh(i) = result(8);

    end
end

plot(SSS, tb_h)

title('Brightness Temperature vs Sea Surface Salinity')
xlabel('SSS (PSU)')
ylabel('T_b (K)')
figure()

plot(SSS, pin_a)
title('Input Power at Detector vs Sea Surface Salinity')
xlabel('SSS (PSU)')
ylabel('P_D (K)')
figure()

plot(SSS, vout_a)
title('Detector Output Voltage vs Sea Surface Salinity')
xlabel('SSS (PSU)')
ylabel('V_{D} (V)')
figure()


plot(pin_a,vout_a)
title('Antenna Voltage vs Input Power Detector')
xlabel('P_d (dBm)')
ylabel('V_{D} (V)')



% plot(ta,ta_n)
% title('Aparent Antenna Temperature (K) vs Antenna Temperature (K) ')
% xlaBMl('Antenna Temperature (K) ')
% ylaBMl('Aparent Antenna Temperature (K)')


function volt_calc = volt_function(T,SSS)


SST = T + 273.15; %Sea Surface Temperature in Kelvin

Tb_h = 0;  %Brightness Temperature (Horizontal Pol)
eia = 0; %Earth Angle of Incidence
volt = 0;
e_inf = 4.9; %permitivitty of free space
e_0 = 8.8412e-12;
freq = 1.413e9; % L-Band frequency

delta = 25 - T;
e_1 = (87.134 - 0.1949*T - 0.01276*T.^2 + 0.0002491*T.^3)*(1 + 1.1613e-5*T*SSS - 0.003656*SSS + 3.21e-5*SSS.^2 - 4.232e-7*SSS.^3);
tau = (1.1109e-10 - 3.824e-12*T + 6.398e-14*T.^2 - 5.096e-16*T.^3)*(1 + 2.282e-5*T*SSS - 7.638e-4*SSS - 7.760e-6*SSS.^2 + 1.105e-8*SSS.^3)/(2*pi);
BMta = 2.033e-2 + 1.266e-4*delta + 2.464e-6*delta.^2 - SSS*(1.849e-5 - 2.551e-7*delta + 2.551e-8*delta.^2);
cond = SSS*(0.18252 - 0.0014619*SSS + 2.093e-5*SSS.^2 - 1.282e-7*SSS.^3)*exp(-delta*BMta); % Conductivty

perm_sw = e_inf + (e_1 - e_inf)/(1 - 2i*pi*freq*tau) + (cond*1i)/(2*pi*freq*e_0); % Permitivity of sea water

gamma_h = (cos(eia) - sqrt(perm_sw - sin(eia).^2))/(cos(eia) + sqrt(perm_sw - sin(eia).^2)); %reflection coefficient

e0_h = 1 - abs(gamma_h).^2; %Smooth Ocean Surface Emissivity (ver pol) %Smooth Ocean Surface Emissivity (hor pol)

e_ocean_h  = e0_h;

Tb_h = e_ocean_h*SST;



NF1 = 10.^(0.35/10); % Noise Figure of LNA 1
NF2 = 10.^(2.29/10); % Noise Figure of RF Switch
NF3 = 10.^(0.49/10); % Noise Figure of BPF 1.7 GHz
NF4 = 10.^(0.35/10); % Noise Figure of LNA 2
NF5 = 10.^(2.5/10); % Noise Figure of BPF 136 MHz

G1_db = 25.454; % Noise Figure of LNA 1 in dB
G2_db = -2.29; % Noise Figure of RF Switch in dB
G3_db = -0.49; % Noise Figure of LNA 2 in dB
G4_db = 25.454; % Noise Figure of LNA 3 in dB
G5_db = -2.5; % Noise Figure of BPF in dB

G1 = 10.^(25.454/10); % Noise Figure of LNA 1
G2 = 10.^(-2.29/10); % Noise Figure of RF Switch
G3 = 10.^(-0.49); % Noise Figure of LNA 2
G4 = 10.^(25.454/10); % Noise Figure of LNA 3
G5 = 10.^(-2.5/10); % Noise Figure of BPF


%NF = NF1 + (NF2 - 1)/G1 + (NF3 - 1)/(G2*G1) + (NF4 - 1)/(G1*G2*G3) + (NF5 - 1)/(G1*G2*G3*G4);

%disp(NF1 + (NF2 - 1)/G1 + (NF3 - 1)/(G2*G1) + (NF4 - 1)/(G1*G2*G3) + (NF5 - 1)/(G1*G2*G3*G4))
disp(G1_db + G2_db + G3_db + G4_db + G5_db)

NF = 10.^(0.5188/10);

%disp(10*log10(NF1 + (NF2 - 1)/G1 + (NF3 - 1)/(G2*G1) + (NF4 - 1)/(G1*G2*G3) + (NF5 - 1)/(G1*G2*G3*G4)))
%disp(G1_db + G2_db + G3_db + G4_db + G5_db)
G_db = 49.1235;
%G_db = 51.5;
G = 10.^(G_db/10);


L = 10.^(0.35/10);
Tp = 290;
Trec = (NF - 1)*Tp;



k = 1.380649e-23;
bw = 136e6;

BM = 0.8808; % eff main beam


p_tbh = 10*log10(1000*(Tb_h + 75.0884+ Trec)*k*bw) + G_db;


L_cable = 10.^(0.35/10); % Perdidas cable

T_cable = (L_cable - 1)*Tp; % temperatura cable

BL = 1; % antenna radiation efficiency


Ta = BL*Tb_h*BM  + Tp*(1 - BL); %ecuacion 4.62 de Ulaby

trec_l = (L_cable - 1)*Tp + L_cable*Trec;

p_in = ((Ta  +  trec_l )*k*bw);
p_in = p_in*1000;

%disp(T_cable)
pin_db = 10*log10(p_in);
pin = pin_db + G_db;

%disp( (10.^(-45.6926/10))*0.001/(k*G*bw))

v_a =  4.33478*(-0.0244*pin+ 0.5651) -6.913;
%v_a =  -0.0244*pin+ 0.5651;

T_ref = 175.8;
p_ref = ((T_ref + trec_l)*k*bw)*1000;
pref_db = 10*log10(p_ref);
pref  = pref_db + G_db;
v_ref = -0.0244*p_ref+ 0.5651;

ta = ((10.^(-(v_a +  2.5855)/0.6971) - 10.^(-(v_ref + 2.5855)/0.6971))*0.001)/(G*k*bw) + 296 ;

ta_n = (ta - (1 - BM)*Tp)/BM;

volt_calc = [pref, v_ref, Tb_h, ta, v_a, pin, ta_n, pin_db];

end