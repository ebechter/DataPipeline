function [velocities,XCfine,masksize,AverageXC] = XCor(spectrum,mask, vel0, vel_step, vel_width,start_order,spec_order,plotflag)


% vel_min = vel0-vel_width;
% vel_max = vel0+vel_width;

N = ceil( (vel_width) / vel_step );  % number of steps 

% velocities2 = [vel_min vel_min + (1:N) * vel_step];


upper = linspace(0,vel_width,N);
lower = -fliplr(upper(2:end));
velocities = [lower upper]+vel0; % symmetric velocities about the central velocity 


delta = linspace(0.07,0.07,39);  % set the mask width - can change per order. 
global colors

for ii = start_order:spec_order
    m_l = mask{ii}(:,1)-delta(ii);
    m_h = mask{ii}(:,1)+delta(ii);
  
    if plotflag ==1
    figure
    set(gcf,'position',[1281 82 2560 420])
    hold on
    box on
    ylim([0 1.2])
    plot(spectrum.wavelength{ii},spectrum.counts{ii},'-k','linewidth',2)
    line([mask{ii}(:,1) mask{ii}(:,1)]', ylim,'Color',[colors{7} 0.8],'LineWidth',2)
%     line([m_h m_h]', ylim)  
    end
    
    masksize(ii)=size(mask{ii}(:,1),1);
        for jj = 1:length(velocities)
            
            XCfine(ii,jj) = CCFcos(m_l,m_h,spectrum.wavelength{ii},spectrum.counts{ii},...
                velocities(jj),mask{ii}(:,2));

%                 velocities(jj),ones(size(mask{ii}(:,1))));

        end
        
end
% norm_masksize = (masksize'-min(masksize))./(max(masksize)-min(masksize));
% weight=repmat(norm_masksize,1,size(XCfine,2));
% WAverageXC = sum(XCfine.*weight,1)./size(XCfine,1);



AverageXC = sum(XCfine,1)./(size(XCfine,1));

% testW = WAverageXC./max(WAverageXC);
% testA = AverageXC./max(AverageXC);


end
