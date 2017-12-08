function make_streamer_input_DYCOMS(zcase)
% Create input files for running Streamer using sounding data for Sc clouds
% Only the DYCOMS case is in the file
%
% zcase determines the resolution in the vertical profile
%    0: Δz=105m
%    1: Δz=30m
%    2: Δz is adaptive, finer in the cloud and above cloud
%
% get Streamer at http://stratus.ssec.wisc.edu/streamer/
% 2017 Monica Zamora, SRAF at UCSD solar.ucsd.edu

%% Sounding data & cloud properties
% From Stevens et al. (2005) http://doi.org/10.1175/MWR2930.1
date_i=datetime(2001,7,15,8,0,0); %UTC time, midnight in local time

z=0:1:3000; % in m
zi=840; ni=find(z==zi,1);
theta_l(z<=zi)=289; theta_l(z>zi)=297+(z(z>zi)-zi).^(1/3); %in K
q_t(z<=zi)=9*1e-3; q_t(z>zi)=1.5*1e-3; %in g/kg

T_sfc=290.4; rho_sfc=1.22;
p_sfc=rho_sfc*287*T_sfc; %in Pa
H=287*T_sfc/9.81; %in m
p=p_sfc*exp(-z/H); rho=rho_sfc*exp(-z/H); %in kg/m3
[T,q_l,q_sat]=get_T_ql_qs(z(1),z,theta_l,q_t,p);


q_v=q_t-q_l;
RH=q_v./q_sat*100;
nb=find(q_l>0,1);
z=z/1e3; zi=zi/1e3; %in km
zb=z(nb); %in km
h=zi-zb; %in km
T_top=T(ni);
LWP=trapz(z,rho.*q_l)*10^3; %LWP not used
LWC=q_l(ni)*rho(ni)*1e3; %LWC at top
PW=-trapz(p,q_t)/1000/9.81*100; %precipitable water in cm
tau=3*LWP/2/1000/1e-5;

%% Reduced system (100 max levels in streamer)
switch zcase
    case 0
        zs=(0:105:3000)/1e3; %00
    case 1
        zs=(0:30:2970)/1e3; %01
    case 2
        ncld=ni-nb; n2=round(ncld/10); %02
        zBL=[z(1),z(round(nb*(1:2)/3)),z(nb-5)];
        zUP=[z(round(ni+ncld+(end-ni-ncld)*(1:5)/5))];
        if n2>=1
            zCLD=z(nb:n2:ni+ncld);
        else
            zCLD=z(nb:ni+ncld);
        end
        zs=[zBL,zCLD,zUP];
end
f=ismember(z,zs);
q_v=q_v(f); p=p(f); T=T(f); RH=RH(f);
z=zs;
%% Create input file
filename_in=['DYCOMS/DYCOMS',num2str(zcase,'%02d'),'.inp'];
fileID = fopen(filename_in,'w');

%% Options section
fprintf(fileID,'OPTIONS \n');
fprintf(fileID,'.TRUE. \n'); %FLUXES compute fluxes?
fprintf(fileID,'.FALSE. \n'); %IR106 use the band 106(3.7um)?
fprintf(fileID,'.TRUE. \n'); %CLDFRC cloud forcing?
fprintf(fileID,'2 2 \n'); %NSTR* SW and LW number of streams
fprintf(fileID,'0 \n'); %NCOEF how many Legendre coeffs
fprintf(fileID,'.TRUE. \n'); %GASABS gaseous absorption?
fprintf(fileID,'.TRUE. \n'); %RAYLISHRT include Rayleigh scattering?
fprintf(fileID,'2 \n'); %ALBTYPE albedo control: model to be scaled by a number
fprintf(fileID,'4 \n'); %EMISSTYPE emissivity control: e for all bands 
fprintf(fileID,'2 .TRUE. \n'); %Profile mid-lat summer, extend to 100km
fprintf(fileID,'4 1 \n'); %Aerosol model: maritime, background tropospheric aerosol
fprintf(fileID,'1 2 3 1 2 2 3\n'); %z(km) T(C) RH(%) zi(km) h(km) λ(Streamer#)
fprintf(fileID,'4 \n'); %Output levels: all
fprintf(fileID,'.TRUE. \n'); %Descriptive output

filename_out=['DYCOMS',num2str(zcase,'%02d'),'.des'];
fprintf(fileID,[filename_out,' \n']); %Output name
fprintf(fileID,'.TRUE. \n'); %User customized output as in writeusr.f
filename_out2=['DYCOMS',num2str(zcase,'%02d'),'.out'];
fprintf(fileID,[filename_out2,' \n']); %
fprintf(fileID,' \n'); % No weights
fprintf(fileID,'.FALSE. \n'); % Won't read cloud optical properties
fprintf(fileID,' \n'); %
fprintf(fileID,' \n'); %

%% Case section
fprintf(fileID,'CASE \n'); %
fprintf(fileID,['DYCOMS II RF01 \n']); %Title
fprintf(fileID,[datestr(date_i,'yy mm dd'),' 8.0 32.85 -117.11 -99 \n']); %Date time and lat-long
fprintf(fileID,' \n'); % Nothing for FLUXES TRUE
fprintf(fileID,'1 129 \n'); %Bands
fprintf(fileID,' 0.05 1 1 1\n'); %Albedo 0.05 w/1 sfc: water
fprintf(fileID,[num2str(T_sfc-273.15),' 1.0 \n']); %Tsrf and emissivity
fprintf(fileID,['1 1 0.9 ',num2str(T_top-273.15),' ',num2str(zi),' ',num2str(h),' ',num2str(tau),' 1 0 10 ',num2str(LWC),' 999 999 999\n']); %Cloud info w/only 1 cloud properties
fprintf(fileID,' \n'); %No cloud overlap
fprintf(fileID,['1 1 1 1 0 0 ',num2str(length(z)),' 1.0 1.0 1.0 1.0 1.0 1.0 \n']); % We read z,p,T,RH profiles
fprintf(fileID,'%.3f %.1f %.1f %.2f \n',[z(end:-1:1);p(end:-1:1)/100;T(end:-1:1)-273.15;RH(end:-1:1)]); %put sounding profiles
fprintf(fileID,[num2str(z(1)),' ',num2str(PW),' \n']); %Sfc height and precipitable water(cm) 

%% Close the file
fclose(fileID);


end
