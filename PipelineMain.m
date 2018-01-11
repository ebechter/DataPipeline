%% Modular Pipeline Main Script  


% clean up before running
clear all
close all

% add all folders below current directory 
addpath(genpath(pwd))

% global custom color list
global colors
scienceframe = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=scienceframe(ii,:);
end
colors{8} =[0.175 0.175 0.175];


% Pipeline options 


Recalibrate = 0; % wavelength calibration = 1, load old solution = 0

RetraceOrders = 1; % find spectral orders again = 1, load old map = 0

ProcessDirectory = 1; % process an entire directory of files = 1, load single file = 0 




% ----------- Image Reduction -----------%

% This section deals with basic image reduction
% Needs dark frame and flat field directories. 
% Not currently used. 

flat = fitsread ('S:\Pipeline\PhotonNoisePipeline\Flat\Format12.22.fits');



% ----------- Wavelength Calibration -----------%


if Recalibrate == 1 

% [wave_cal] = WavelengthCalibration(wavecalframe,etalonlinelist);

else
    
% This skips re-calibration and loads in coefficients used in the simulator
coeffs = 'polycoeffs12.22.mat';
[wave_cal] = ReloadWavelengthSolution(coeffs);
end 


% ----------- Spectral Extraction -----------%

% Load appropriate files


if ProcessDirectory == 1
    
    % Run a directory
    sciencepath = 'S:\Pipeline\ModularPipeline\DetectorFrames\';
    file_id = strcat(sciencepath,'*.fits');
    allFiles = dir(file_id);
    fname = {allFiles.name};
    
else
    
    sciencepath = 'S:\Pipeline\ModularPipeline\DetectorFrames\';
    fname = {'.fits'};
    
end


if RetraceOrders == 1
    
    % retrace orders
    trace = BasicOrderIdentification(flat);
    save tracefile1.mat trace
    
else
    
    tracefile = '.mat';
    load(tracefile);
    
end

    window =5;
    norm=0;
    
    
    
    for ii = 1:size(fname,2)  % Work on all files

    
    scienceframe = fitsread([sciencepath fname{ii}]);
    
    
    % Extract all orders here
    ex_orders{ii} = BasicOrderExtraction(scienceframe,trace,window,norm);
    end
    
%     ex_science = HOrderExtraction(science,window,0); % for physically
%     flattened orders only

return





% ----------- Mask Building  -----------%

load mask16
mask(:,1)=mask(:,1)*10;
offset = 0; 





% trace = BasicOrderIdentification(flat);
orderflat = HOrderExtraction(flat,window,0);

% trace = trace(:,tracenum:3:end);
col = 4;
sciencepath = 'S:\Pipeline\PhotonNoisePipeline\Polarization\100\';
file_id = strcat(sciencepath,'*.fits');
allFiles = dir(file_id);
fname = {allFiles.name};
% path = '';
% fname = {'S:\Pipeline\PhotonNoisePipeline\Polarization\1\Polarized1pfrac50.fits'};
start_order = 1;
spec_order =  36;

for ii = 1:size(fname,2)  % Work on all files
    clear msk
    clear sciencespec
    clear science
    
    fprintf('\nWorking on file %s\n',fname{ii})
    

    xx = 0.5:4095.5; % make x vector for order polynomial solution
    
    % Build spectrum and mask for each order... in Angstroms
    
    for jj =  start_order:spec_order
        %         (polyval(wave_solution(jj,:),xx)./1e3)'
        sciencespec{jj} = [zem(:,jj),scienceframe(:,jj)./(orderflat(:,jj)./max(orderflat(:,jj)))];
        
        sciencespec{jj} = LinDetrend(sciencespec{jj});
        
        % Can doppler shift spectrum right before cross correlation (for testing)
        RV = [0 0]*1e3;
        c = 2.9979245800 * 10^8;  % Speed of light [m/s] according to NIST - http://physics.nist.gov/cgi-bin/cuu/Value?c
        beta = RV ./ c;
        delta = sqrt((1 + beta) ./ (1 - beta));
        lower_bound = sciencespec{jj}(1,1)* delta(2);
        upper_bound = sciencespec{jj}(end,1)* delta(1);
        
        ind =  mask(:,1) >= lower_bound  &  mask(:,1) <= upper_bound;
        msk{jj} = mask(ind,:);
        
    end
    
    msk{1}([16],:) = [];
    msk{2}([2 3 4],:) = [];
    msk{3}([4 5 8],:) = [];
    msk{4}([6 7 11 12 22 23 26 27 38 39 44 45 60 61 64 65 66 70 71],:) = [];
    msk{5}([1 2 3 8 9 16:19 20 21 24 25 26 27 33 46 50 55 56 61],:) = [];
    msk{6}([6 7 12 13 25 26 32 33 43 44 45],:) = [];
    msk{7}([1:7 15 16 22 23 24 46 47 48 49 50],:) = [];
    msk{8}([12 13 14 25 26 27 28 30 31 35 36],:) = [];
    msk{9}([8 12 13 14 17 19 20 21 25 26 27 28 29 31 34 45 46 58 59 62 63 64 65 66 74 75 78 79 ],:) = [];
    msk{10}([3,7,8,9,10,12,13,14,15,19,20,26,27,28,29,34,44,45,46],:) = [];
    msk{11}([5 7 9 13 14 19:20 27 31 33 36 37 39 42 44],:) = [];
    msk{12}([5,25,26,29,34,35,36,37,39,40,42],:) = [];
    msk{13}([1 4 5 19:22 25 26 32 36 37 42 43 46 47],:) = [];
    msk{14}([4 8 10 15 18:20 22 25:28 35 37],:) = [];
    msk{15}([5 6 12 14 15 17],:) = [];
    msk{16}([11 22 24 25 31 32 33 ],:) = [];
    msk{17}([1 19 25 26],:) = [];
    msk{18}([1 15 16 25],:) = [];
    msk{19}([2 5 12 13],:) = [];
    msk{20}([11 15 16],:) = [];
    msk{21}([5 19],:) = [];
    msk{22}([3 4 9 10 14 15],:) = [];
    msk{23}([10 13 15],:) = [];
    msk{24}([],:) = [];
    msk{25}([5 6 9],:) = [];
    msk{26}([3 5 6],:) = [];
    msk{27}([1 5 3],:) = [];
    msk{28}([1 2],:) = [];
    msk{29}([],:) = [];
    msk{30}([],:) = [];
    msk{31}([],:) = [];
    msk{32}([],:) = [];
    msk{33}([],:) = [];
    msk{34}([1],:) = [];
    msk{35}([1 2],:) = [];
    msk{36}([],:) = [];
    
    % Set initial search parameters for XCor.. km/s
    vel0 = 0;
    vel_step = 0.4;
    vel_width = 35;
    
    [XCvelCoarse, XCC, masksize, AVE] = XCor(sciencespec,msk, vel0, vel_step, vel_width,1,36);
    
    [~,temp] = min(AVE);
    vel0 = XCvelCoarse(temp);
    
    [RoughRV] = FitOrderXCor(1,1,XCvelCoarse,vel0,AVE,masksize,0);
    
    vel0 = RoughRV;  % Refine
    vel_step = 0.1;
    vel_width = 25;
    
    % Cross Correlate from start_order to spec_order
    [XCvelocities, XCF, masksize, AV_XC] = XCor(sciencespec,msk, vel0, vel_step, vel_width,start_order,spec_order);
    
    % Fit order XC functions individually
    [order_RV(ii,:)] = FitOrderXCor(start_order,spec_order,XCvelocities,vel0,XCF,masksize,0);
    
    % Fit the average XC function
    gstartpoints = [200 vel0 1 0];
    [qual,aveans(ii)] = GaussianFit(XCvelocities,AV_XC',gstartpoints,1,colors{rem(ii,length(colors)-1)+1},3);
    fprintf('Average Fit RV residual (m/s): %.5f, R^2: %.6f\n',aveans(ii)*1000-offset,qual)
    
    % Plot order performance
    if OrderPerformance ==1
        figure(1000)
        hold on
        h(ii)=plot(order_RV(ii,:)*1000*100-offset*100,'-','color',colors{rem(ii,length(colors)-1)+1},'marker','.','markersize',14);
        box on
        grid on
        xlabel('order number','fontsize',14)
        ylabel('Velocity residual (cm/s)','fontsize',14)
    end
    [HeaderCell] = fitshead([sciencepath fname{ii}]);
    idx = find(strcmp([HeaderCell],'SMFPOSX')); % single line engine
    smfx(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'SMFPOSY')); % single line engine
    smfy(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'SMFPOSZ')); % single line engine
    smfz(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'ADCPOS')); % single line engine
    adc(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'RHO')); % single line engine
    Rho(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'AIRMASS')); % single line engine
    Airmass(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'DOP')); % single line engine
    dop(ii) = HeaderCell{idx,2};
    idx = find(strcmp([HeaderCell],'PFRAC')); % single line engine
    pfrac(ii) = HeaderCell{idx,2};
end
residualcm = 100*(aveans*1000-offset);

figure(1)
hold on
Pol = [pfrac',dop',residualcm'];
Pol = sortrows(Pol,1);
% plot(Pol(:,1),Pol(:,3),'.','color',colors{1})
plot(Pol(:,1),Pol(:,3),'linewidth',2,'color',colors{col})
grid on
box on
ax = gca;
ax.LineWidth = 1.5;

%Fiber alignment
radialpos = sqrt(smfx.^2+smfy.^2);

%Polarization

dop;
pfrac;
return

% AO and ADC
figure
adcrv(:,1) = Airmass;
adcrv(:,2) = residualcm;
adcrv(:,3) = Rho;
test = sortrows(adcrv,1);
plot(180*asec(test(:,1))/pi,test(:,3),'-','linewidth',1.5)
grid on
box on
ax = gca;
ax.LineWidth = 1.5;
ylabel('RV offset (cm/s)')
xlabel('Zenith distance (degrees)')

figure
adcrv(:,1) = adc;
adcrv(:,2) = residualcm;
adcrv(:,3) = Rho;
test = sortrows(adcrv,1);
plot([0:5:60],test(:,2),'-','linewidth',1.5)
grid on
box on
ax = gca;
ax.LineWidth = 1.5;
ylabel('RV offset (cm/s)')
xlabel('Zenith distance (degrees)')


