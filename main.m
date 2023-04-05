clc
clear all
% close all

%this code considers the size distribution of spherical particles

%solution of radiative transfer equation in one layer pigmented plane parallel medium
%ray is incident from air to coating. coating is coated on a substrate
%substrate could be air or other material such as silver, glass etc.
%the code estimates spectral hemispherical reflectance, transmittance and absorptance
%can handle; independent scattering, boundary reflections, absorption in
%medium. can't handle; coherent backscattering, dependent scattering and polarized ray tracing.
%while calculating the scattering direction the code uses cumulative inverse
%relation of exact single scattering pahse function or henyey greenstein 
%function approximation depending on choice 

lamda=(400:100:20000)*10^-9; %freespace wavelength of incident ray in meter 
thickness=1*10^-3;  %thickness of coating in meter 
lamda_um = lamda*10^6;

f_v=0.04; %volume fraction of air bubbles

polar_angle=linspace(0,89.99999,9); %incident angles. 0 = perpendicular to slab face. 90 parallel and should be avoided.


r_avg=5;%um
std_dev=0.01 ;
r_vector_um=(1:0.2:15); %size range to be considered
weight_vector = lognpdf(r_vector_um,log(r_avg),std_dev); %the distribution
r_vector=r_vector_um*10^-6;
weight_vector=weight_vector/trapz(r_vector,weight_vector); 

figure %show the size distribution
plot(10^6*r_vector,weight_vector,'LineWidth',2)
ylabel('Frequency')
xlabel('Radius [\mum]')

[n_medium, k_medium] = PDMS_nk(lamda*10^6); %Compute n, k of PDMS at lamda in um

n_pigment=ones(length(lamda),1); %real refractive index of substrate
k_pigment=zeros(length(lamda),1); %imaginary refractive index of substrate
n_substrat=ones(length(lamda),1); %real refractive index of substrate
k_substrat=zeros(length(lamda),1); %imaginary refractive index of substrate

photon_number=10^4; %number of rays that will be traced, higher the number more accurate the result
n_cdf_random=1000; %how many pieces will be between (0,1) in random number relation with inverse cumulative function, higher the number accurate the phase function. useless if HG is used
nang_gid=1000; %how many pieces will be between (0,pi) in cumulative distribution function, higher the number more accurate the result. useless if HG is used


lamda_nm=lamda*10^9; %for plot
polar_angle_rad=polar_angle*pi/180;


%initialize arrays
ref_lamda=zeros(length(lamda),length(polar_angle));
tra_lamda=zeros(length(lamda),length(polar_angle));
abs_lamda=zeros(length(lamda),length(polar_angle));
C_sca_sigma=zeros(length(lamda),length(r_vector));
C_abs_sigma=zeros(length(lamda),length(r_vector));
g_arr_sigma=zeros(length(lamda),length(r_vector));

% calculate surface reflection from air to medium. 
% medium to air and medium to substrate is calculated within snell.m.
% air to medium is calculated here seperately since we need refraction angle
teta_prime=zeros(length(lamda),length(polar_angle));
sur_reflection=zeros(length(lamda),length(polar_angle));
for i=1:length(lamda)
    for j=1:length(polar_angle_rad)
        teta_prime(i,j)=F_fresnel_2(n_medium(i),k_medium(i),polar_angle_rad(j))*180/pi;
        cos_teta=cosd(polar_angle(j));
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=1/(n_medium(i)-1i*k_medium(i));
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        carpan1=cos_teta/cos_x2;
        carpan3=cos_x2/cos_teta;
        E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=(carpan3-carpan2)/(carpan3+carpan2);
        R_orth=E_orth*conj(E_orth);
        reflectance=real(R_parallel+R_orth)*0.5;
        sur_reflection(i,j)=reflectance;
    end
end

for z=1:length(r_vector)
    Area=pi*r_vector(z)^2;
    for i=1:length(lamda)
        x=2*pi*r_vector(z)*n_medium(i)/lamda(i);
        m=(n_pigment(i)+1i*k_pigment(i))/n_medium(i);
        fonksiyon=Mie(m,x);
        Qabs=fonksiyon(3);
        Qsca=fonksiyon(2);
        C_sca_sigma(i,z)=Qsca*Area;
        C_abs_sigma(i,z)=Qabs*Area;
        g_arr_sigma(i,z)=fonksiyon(5);
    end
end

r_3=trapz(r_vector,(weight_vector.*r_vector.^3)');
V_avg=(4/3)*pi*r_3;
C_sca=trapz(r_vector,(weight_vector.*C_sca_sigma)');
C_abs=trapz(r_vector,(weight_vector.*C_abs_sigma)');

sigma_s = C_sca*f_v/V_avg;
kappa = C_abs*f_v/V_avg + (1-f_v)*4 * pi * k_medium./lamda;
ext_tot = sigma_s + kappa;
g=trapz(r_vector,(weight_vector.*g_arr_sigma.*C_sca_sigma)')./trapz(r_vector,(weight_vector.*C_sca_sigma)');
scat_prob=sigma_s./ext_tot;

tic
%loop the monte carlo code for lamda and polar_angle
for j=1:length(polar_angle)
    parfor i=1:length(lamda)
        [ref_lamda(i,j),tra_lamda(i,j),abs_lamda(i,j)]=monte_carlo(photon_number,sur_reflection(i,j),cosd(teta_prime(i,j)),thickness,scat_prob(i),ext_tot(i),n_medium(i),k_medium(i),n_substrat(i),k_substrat(i),g(i));
    end
end
toc
figure %draw normal to diffuse R, T and A for normal incidence (first index in my case)
plot(lamda_um,ref_lamda(:,1),'--r',lamda_um,tra_lamda(:,1),'--k','LineWidth',2)
ylim([0 1])
xlim([0 20])
legend('R_n_h','T_n_h','Location', 'Best')
xlabel('Wavelength [\mum]')
ylabel({'Normal to Hemispherical';'Reflectance, Transmittance, Absorptance'})
error('end')
T_sur=300;
T_amb=300;
BB_Tsur_l=I_bb(lamda,T_sur)';
emittance_cal=abs_lamda(:,1);


trans_atm_calc=emis_atm_new(lamda)';



figure('Renderer', 'painters', 'Position', [500 300 428 420]) % starting point and height - width of the frame
hAx=gca;
area(lamda_um,trans_atm_calc,'LineStyle','none')
hold on
plot(lamda_um,BB_Tsur_l/max(BB_Tsur_l),'-r',lamda_um,emittance_cal,'--k','LineWidth',2)
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
ylim([0 1])
xlim([2 20])
hLg=legend('Atmospheric transmissivity','Blackbody intensity at 300K','Calculation','Experiment_s','Experiment_r','Location', 'southeast')
hLg.LineWidth=1.5;
hLg.EdgeColor = [0 0 0];
xlabel('Wavelength [\mum]')
ylabel('Spectral normal emittance')
set(gca,'FontSize',13)
set(gca,'XMinorTick','on','YMinorTick','on')
box on
saveas(gcf,'abs.png')
saveas(gcf,'abs.emf')
saveas(gcf,'abs.fig')




BB_Tsur_l=I_bb(lamda,T_sur)';

emis_calc=trapz(lamda,emittance_cal.*BB_Tsur_l)/trapz(lamda,BB_Tsur_l)

emis_calc_window=trapz(lamda,emittance_cal.*BB_Tsur_l.*trans_atm_calc)/trapz(lamda,BB_Tsur_l.*trans_atm_calc)


emittance=1-ref_lamda-tra_lamda;

polar_angle_rad=polar_angle*pi/180;
Solar_l=I_solar(lamda)';

trans_atm_pure=emis_atm_new(lamda)';
emit_atm=1-trans_atm_pure.^(1./cos(polar_angle_rad));

BB_Tamb_l=I_bb(lamda,T_amb)';
BB_Tsur_l=I_bb(lamda,T_sur)';
emittance_rad=trapz(lamda,BB_Tsur_l.*emittance,1);
emittance_atm=trapz(lamda,BB_Tamb_l.*emittance.*emit_atm,1);

P_rad=trapz(polar_angle_rad,2*cos(polar_angle_rad).*sin(polar_angle_rad).*emittance_rad)
P_atm=trapz(polar_angle_rad,2*cos(polar_angle_rad).*sin(polar_angle_rad).*emittance_atm)
P_sol=trapz(lamda,Solar_l.*emittance(:,1))
P_net=P_rad-P_atm- P_sol


