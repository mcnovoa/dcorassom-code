
NF1 = 10.^(1.7/10); % Noise Figure of LNA 1
NF2 = 10.^(3.09/10); % Noise Figure of RF Switch
NF3 = 10.^(1.7/10); % Noise Figure of BPF 1.7 GHz
NF4 = 10.^(1.7/10); % Noise Figure of LNA 2
NF5 = 10.^(4.25/10); % Noise Figure of BPF 136 MHz
NF6 = 10.^(6.5/10); % Noise Figure of BPF 136 MHz

G1_db = 18.606; % Noise Figure of LNA 1 in dB
G2_db = -3.099; % Noise Figure of RF Switch in dB
G3_db = 18.505; % Noise Figure of LNA 2 in dB
G4_db = 18.606; % Noise Figure of LNA 3 in dB
G5_db = -4.25; % Noise Figure of BPF in dB
G6_db = 26.48; % Noise Figure of BPF in dB

G1 = 10.^(G1_db/10); % Noise Figure of LNA 1
G2 = 10.^(G2_db/10); % Noise Figure of RF Switch
G3 = 10.^(G3_db/10); % Noise Figure of LNA 2
G4 = 10.^(G4_db/10); % Noise Figure of LNA 3
G5 = 10.^(G5_db/10); % Noise Figure of BPF
G6 = 10.^(G6_db/10);

NF_db = NF1 + (NF2 - 1)/G1 + (NF3 - 1)/(G2*G1) + (NF4 - 1)/(G1*G2*G3) + (NF5 - 1)/(G1*G2*G3*G4) + (NF6 - 1)/(G1*G2*G3*G4*G5);

NF = 10*log10(NF_db);

k = 1.380649e-23;
bw = 8e9;
G_sys= G1_db + G2_db + G3_db + G4_db + G5_db + G6_db;

%disp(0.02744*(3.5115) + 0.9788)

L = 1;
Tp = 290;
Trec = (NF - 1)*Tp;

L_cable = 10.^(1.1/10);


trec_l = (L_cable - 1)*Tp + L_cable*Trec;

disp(NF)
disp(G_sys)
disp(trec_l)