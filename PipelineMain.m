%% Modular Pipeline Main Script
% clean up before running
clear all
clc
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

RetraceOrders = 0; % find spectral orders again = 1, load old map = 0

ProcessDirectory = 1; % process an entire directory of files = 1, load single file = 0

trace = 1;
start_order = 1;
spec_order = 36;
load mask16.mat
mask(:,1)=mask(:,1)*10;
offset = 0;
% ----------- Image Reduction -----------%

% This section deals with basic image reduction
% Needs dark frame and flat field directories.
% Not currently used.

flat = fitsread ('/Volumes/Software/Pipeline/PhotonNoisePipeline/Flat/Format12.22.fits');



% ----------- Wavelength Calibration -----------%


if Recalibrate == 1
    
    % [wave_cal] = WavelengthCalibration(wavecalframe,etalonlinelist);
    
else
    
    % This skips re-calibration and loads in coefficients used in the simulator
    coeffs = 'polycoeffs12.22.mat';
    [wave_cal] = ReloadWavelengthSolution(coeffs,trace);
end


% ----------- Spectral Extraction -----------%

% Load appropriate files


if ProcessDirectory == 1
    
    % Run a directory
    sciencepath = '/Volumes/Software/Pipeline/ModularPipeline/DetectorFrames/Coma/';
    file_id = strcat(sciencepath,'*.fits');
    allFiles = dir(file_id);
    fname = {allFiles.name};
    
else
    
    sciencepath = '/Volumes/Software/Pipeline/ModularPipeline/DetectorFrames/';
    fname = {'ScienceTraceTestFull.fits'};
    
end


if RetraceOrders == 1
    
    % retrace orders
    trace = BasicOrderIdentification(flat);
    save tracefile1.mat trace
    
else
    
    tracefile = 'tracefile1.mat';
    load(tracefile,'-mat');
    trace = trace(:,1:3:end);
end

window =15;
norm=0;



for ii = 1:size(fname,2)  % Work on all files
    
    
    scienceframe = fitsread([sciencepath fname{ii}]);
%     scienceframe = load([sciencepath fname{ii}],'-mat');
    
    % Extract all orders here
    ex_orders{ii} = BasicOrderExtraction(scienceframe,trace,window,norm);
    
    
    %     ex_science = HOrderExtraction(science,window,0); % for physically
    %     flattened orders only
    
    xx = 0.5:4095.5; % make x vector for order polynomial solution
    
    
    % Build spectrum and mask for each order... in Angstroms
    
    for jj =  start_order:spec_order
        %         (polyval(wave_solution(jj,:),xx)./1e3)'
        %         figure(1)
        %         hold on
        %
        dlam = diff(wave_cal(:,jj));
        dlam = [dlam(1); dlam];
        
        spectrum.wavelength{jj} = wave_cal(:,jj);
        spectrum.counts{jj} = ex_orders{ii}(:,jj);
        
        %
        %         if jj ==1
        %             normfactor = max(sciencespec{jj}(:,2)./dlam);
        %         end
        
        %./orderflat(:,jj)];
        %         plot(sciencespec{jj}(:,1),(sciencespec{jj}(:,2)./dlam))        % Can doppler shift spectrum right before cross correlation (for testing)
        
        %         sciencespec{jj} = LinDetrend(sciencespec{jj});
        %         plot(sciencespec{jj}(:,1),sciencespec{jj}(:,2),'color',[0 0 0 0.8])        % Can doppler shift spectrum right before cross correlation (for testing)
        RV = [0 0]*1e3;
        c = 2.9979245800 * 10^8;  % Speed of light [m/s] according to NIST - http://physics.nist.gov/cgi-bin/cuu/Value?c
        beta = RV ./ c;
        delta = sqrt((1 + beta) ./ (1 - beta));
        lower_bound = spectrum.wavelength{jj}(1)* delta(2);
        upper_bound = spectrum.wavelength{jj}(end)* delta(1);
        
        %         sciencespec{jj}(:,1) = sciencespec{jj}(:,1).*delta;
        
        ind =  mask(:,1) >= lower_bound  &  mask(:,1) <= upper_bound;
        
        msk{jj} = mask(ind,:);
        %         msk{jj}(:,1) = lam_centers{jj}; % build mask for specific order jj
        
        %         msk{jj}=msk{jj}(SkipPeakNum+1:end-SkipPeakNum); % trim mask by a few peaks
        
    end
    %     maskcopy = msk;
    %     msk{1}([16],:) = [];
    %     msk{2}([2 3 4],:) = [];
    %     msk{3}([4 5 8],:) = [];
    %     msk{4}([6 7 11 12 22 23 26 27 38 39 44 45 60 61 64 65 66 70 71],:) = [];
    %     msk{5}([1 2 3 8 9 16:19 20 21 24 25 26 27 33 46 50 55 56 61],:) = [];
    %     msk{6}([6 7 12 13 25 26 32 33 43 44 45],:) = [];
    %     msk{7}([1:7 15 16 22 23 24 46 47 48 49 50],:) = [];
    %     msk{8}([12 13 14 25 26 27 28 30 31 35 36],:) = [];
    %     msk{9}([8 12 13 14 17 19 20 21 25 26 27 28 29 31 34 45 46 58 59 62 63 64 65 66 74 75 78 79 ],:) = [];
    %     msk{10}([3,7,8,9,10,12,13,14,15,19,20,26,27,28,29,34,44,45,46],:) = [];
    %     msk{11}([5 7 9 13 14 19:20 27 31 33 36 37 39 42 44],:) = [];
    %     msk{12}([5,25,26,29,34,35,36,37,39,40,42],:) = [];
    %     msk{13}([1 4 5 19:22 25 26 32 36 37 42 43 46 47],:) = [];
    %     msk{14}([4 8 10 15 18:20 22 25:28 35 37],:) = [];
    %     msk{15}([5 6 12 14 15 17],:) = [];
    %     msk{16}([11 22 24 25 31 32 33 ],:) = [];
    %     msk{17}([1 19 25 26],:) = [];
    %     msk{18}([1 15 16 25],:) = [];
    %     msk{19}([2 5 12 13],:) = [];
    %     msk{20}([11 15 16],:) = [];
    %     msk{21}([5 19],:) = [];
    %     msk{22}([3 4 9 10 14 15],:) = [];
    %     msk{23}([10 13 15],:) = [];
    %     msk{24}([],:) = [];
    %     msk{25}([5 6 9],:) = [];
    %     msk{26}([3 5 6],:) = [];
    %     msk{27}([1 5 3],:) = [];
    %     msk{28}([1 2],:) = [];
    %     msk{29}([],:) = [];
    %     msk{30}([],:) = [];
    %     msk{31}([],:) = [];
    %     msk{32}([],:) = [];
    %     msk{33}([],:) = [];
    %     msk{34}([1],:) = [];
    %     msk{35}([1 2],:) = [];
    %     msk{36}([],:) = [];
    
    %         msk{2}([2,8,10,12,end],:) = [];
    %         msk{3}([5],:) = [];
    % Set initial search parameters for XCor.. km/s
    %     msk{33}([7 8 9 14 17 18 20 21],:) = [];
    
    %     for mm = 1:size(msk{33},1)
    %         msk{33} = mask(ind,:);
    %
    %         msk{33} = msk{33}(mm,:);
    %
    injected = 0;
    vel0 = 0;
    vel_step = 0.4;
    vel_width = 50;
    
    [XCvelCoarse, XCC, masksize, AVE] = XCor(spectrum,msk, vel0, vel_step, vel_width,start_order,spec_order,0);
    
    [~,temp] = min(AVE);
    vel0 = XCvelCoarse(temp);
    
    [RoughRV] = FitOrderXCor(1,1,XCvelCoarse,vel0,AVE,masksize,0,injected(end));
    
    vel0 = RoughRV;  % Refine
    vel_step = 0.1;
    vel_width = 35;
    
    
    % Cross Correlate from start_order to spec_order
    [XCvelocities, XCF, masksize, AV_XC] = XCor(spectrum,msk, vel0, vel_step, vel_width,start_order,spec_order,0);
    
    
    % Fit order XC functions individually
    [order_RV(ii,:)] = FitOrderXCor(start_order,spec_order,XCvelocities,vel0,XCF,masksize,0,injected(end));
    
    
    % Fit the average XC function
    gstartpoints = [200 vel0 1 0];
    [qual,aveans(ii)] = GaussianFit(XCvelocities,AV_XC',gstartpoints,0,colors{rem(ii,length(colors)-1)+1},3);
    fprintf('Average Fit RV residual (m/s): %.5f, R^2: %.6f\n',aveans(ii)*1000,qual)
    
    
    
    %     fprintf('Average Fit RV residual (m/s): %.5f, R^2: %.6f\n',aveans(ii)*1000-injected(end)-offset,qual)
end

plotvals = sort(aveans)*1000;
figure(1)
hold on
plot(plotvals,'.','markersize',15)

