% This program is used to determine the stress induced between welds on 
% interconnects that results from an increase in temperature
% By Katherine Han
% SunPower
% 9-5-14

% governing equation: delta_T = sigma/(E*delta_alpha) + C *
% sigma_yield/(E*delta_alpha)*(sigma/sigma_yield)^m
% where delta_T is the increase in temperature from the unstrained state 
% in deg C, sigma is the stress in MPa, E is the young's modulus in MPa,
% delta_alpha is the difference in the coefficient of thermal expansions
% between the metal of interest and silicon in K^-1, C is an emperical
% fitting factor, sigma_yield is the yield stress of the material, and m is
% another fitting factor

clear all
close all

% Geometric inputs for interdigitated fingers
t = 0.37:0.1:0.47; % thicknes of fingers (foil thickness) in mm
L = 0.5:0.5:2; % distance between welds on fingers in mm

% Give a range of the t/l ratio to sweep
tl = 0.01:0.01:.05;
tl(6) = 0.0185;
tl(7) = 0.037;
tl(8) = 0.043529;
tl(9) = 0.074;
tl(10) = 0.003803;


% Set up desired temperature range to graph later (above room temperature), where 0 is 25
% degrees C and 75 is 100 degrees C, and desired range of sigmas to sweep
T = 0:175; % increase in T from 25C
sigma = 0:250; %sigma to sweep in MPa

% number of materials in this simulation
num_metals = 2;

% Input material properties for metals, Ramberg-Osgood parameters
        E_Al_1145_0 = 69000; %Young's modulus in MPa
        C_Al_1145_0 = 3.529; % empirical fitting parameter
        m_Al_1145_0 = 8.074; % exponential empirical fitting parameter
        sigma_yield_Al_1145_0 = 43.5199; %yield strenght of the material in MPa
        delta_alpha_Al_1145_0 = 0.0000213; % difference in coefficient of thermal expansion of metal to silicon
        
        E_Cu = 110000; %Young's modulus in MPa
        C_Cu = 4.2248; % empirical fitting parameter
%         m_Cu = 5.6746; % exponential empirical fitting parameter
        m_Cu = 8; % exponential empirical fitting parameter
        sigma_yield_Cu = 126.487; %yield strenght of the material in MPa
        delta_alpha_Cu = 0.0000142; % difference in coefficient of thermal expansion of metal to silicon
    
% Initialize delta T output vector with a row for each type of metal
delT_metal = zeros(num_metals,length(sigma));

% Solve for delta T based on a range of stresses
for sig = 1:length(sigma)
    delT_metal(1,sig) = sigma(sig)/E_Al_1145_0/delta_alpha_Al_1145_0 +...
        C_Al_1145_0*sigma_yield_Al_1145_0/E_Al_1145_0/delta_alpha_Al_1145_0*...
        (sigma(sig)/sigma_yield_Al_1145_0)^m_Al_1145_0;
    delT_metal(2,sig) = sigma(sig)/E_Cu/delta_alpha_Cu +...
        C_Cu*sigma_yield_Cu/E_Cu/delta_alpha_Cu*(sigma(sig)/sigma_yield_Cu)^m_Cu;

end

delT_tl = zeros(num_metals,length(tl));
sigma_c = zeros(length(num_metals),length(tl));
for tl1 = 1:length(tl)
    
    delT_tl(1,tl1) = pi^2*tl(tl1)^2/3/delta_alpha_Al_1145_0;
    sigma_c(1,tl1) = spline(delT_metal(1,:),sigma(:),delT_tl(1,tl1));
    
    delT_tl(2,tl1) = pi^2*tl(tl1)^2/3/delta_alpha_Cu;
    sigma_c(2,tl1) = spline(delT_metal(2,:),sigma(:),delT_tl(2,tl1));

end


% Measure the deflection with changing temperature for all thickness/Length
% combinations, assume a thickness of 

deflection37 = zeros(length(T),length(tl));
deflection135 = zeros(length(T),length(tl));


Length37 = 0.037./tl;
Length135 = 0.135./tl;

delta_alpha_Al_1145_0;

for Len = 1:length(tl)
    for temperature = 1:length(T)
        deflection37(temperature,Len) = 2*Length37(Len)/pi*...
            (delta_alpha_Al_1145_0*T(temperature)+(delta_alpha_Al_1145_0*T(temperature))^2/2)^0.5;
        deflection135(temperature,Len) = 2*Length135(Len)/pi*...
            (delta_alpha_Al_1145_0*T(temperature)+(delta_alpha_Al_1145_0*T(temperature))^2/2)^0.5;
    end
end
        



figure
plot(delT_metal(1,:),sigma,'b')
hold on
plot(delT_metal(2,:),sigma,'g')
hold on
plot(delT_tl(1,1:5),sigma_c(1,1:5),'ob','MarkerSize',10)
hold on
plot(delT_tl(2,1:5),sigma_c(2,1:5),'xg','MarkerSize',10)
hold on
plot(delT_tl(1,6),sigma_c(1,6),'+r','MarkerSize',10)
hold on
plot(delT_tl(2,6),sigma_c(2,6),'+k','MarkerSize',10)
hold on
plot(delT_tl(1,7),sigma_c(1,7),'*r','MarkerSize',10)
hold on
plot(delT_tl(2,7),sigma_c(2,7),'*k','MarkerSize',10)
hold on
plot(delT_tl(1,8),sigma_c(1,8),'sr','MarkerSize',10)
hold on
plot(delT_tl(2,8),sigma_c(2,8),'sk','MarkerSize',10)
hold on
plot(delT_tl(1,9),sigma_c(1,9),'dr','MarkerSize',10)
hold on
plot(delT_tl(2,9),sigma_c(2,9),'dk','MarkerSize',10)
hold on
plot(delT_tl(1,10),sigma_c(1,10),'vr','MarkerSize',10)
hold on
plot(delT_tl(2,10),sigma_c(2,10),'vk','MarkerSize',10)
legend('Aluminum 1145-O','Copper','Critical Stress, Al','Critical Stress, Cu','Al L=2 mm, t=0.037mm','Cu L=2 mm, t=0.037mm'...
    ,'Al L=1 mm, t=0.037mm','Cu L=1 mm, t=0.037mm','Al L=0.85 mm, t=0.037mm','Cu L=0.85 mm, t=0.037mm'...
    ,'Al L=0.5 mm, t=0.037mm','Cu L=0.5 mm, t=0.037mm','Al L=35.5 mm, t=0.135mm','Cu L=2 mm, t=0.135mm','FontSize',14,'Location','best')
xlim([0 400])
xlabel('Temperature (C)','FontSize',16)
ylabel('Stress (MPa)','FontSize',16)
% legend('Aluminum','Copper')


% for k = 1:length(Length37)
%     
%     Length37str(k) = (num2str(Length37(k)));
%     Length135str(k) = num2str(Length135(k));
% end
% 

figure
plot(T,deflection37(:,1),'-r',T,deflection37(:,2),'-m',T,deflection37(:,3),'-g',...
    T,deflection37(:,4),'-c',T,deflection37(:,5),'-b',T,deflection37(:,6),'-k',T,deflection37(:,7),'-b',...
    T,deflection37(:,8),'-b',T,deflection37(:,9),'-b',T,deflection37(:,10),'-b')
title('Deflection for 37 um foil for weld distances in mm')
xlabel('Temperature (C)')
ylabel('Deflection (mm)')
xlim([0 175])
ylim([0 0.2])
% legend(Length37str')


figure
plot(T,deflection135(:,1),':r',T,deflection135(:,2),':m',T,deflection135(:,3),':g',...
    T,deflection135(:,4),':c',T,deflection135(:,5),':b',T,deflection135(:,6),':k',T,deflection135(:,7),':k',...
    T,deflection135(:,8),':k',T,deflection135(:,9),':k',T,deflection135(:,10),':k')
xlim([0 175])
ylim([0 0.5])
title('Deflection for 135 um foil for weld distances in mm')
xlabel('Temperature (C)')
ylabel('Deflection (mm)')
% legend(Length135str')


%  for thick = 1:length(t)
%     for weld_distance = 1:length(L)
%         for sigma_sweep = 1:length(sigma)
%                 metal = met;
%                 
       
% set up cases for each type of metal finger 

