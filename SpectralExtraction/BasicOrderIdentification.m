function [trace,isOrderComplete] = BasicOrderIdentification(flat)

%% Quick look 
% figure(1000)
% imagesc(flat)
global colors
% Define extraction window in pixels and number of samples 
lowerlim = 1; % start at pixel 500
upperlim = 4096; % end at pixel 3200 
samples = 50; % 20 samples 

xx=linspace(lowerlim,upperlim,samples); 
peak = zeros(36*3,length(xx));
loc = zeros(36*3,length(xx));

isOrderComplete = ones(1,36);

% Smooth crosscut using a gaussian, record peaks and locations  
for ii = 1:samples
  
    b = smoothts(flat(:,floor(xx(ii))),'g',5,0.65);
    [pk,lk] = findpeaks(b);
%     if length(pk) ~= 117
%         
%         n = 117 - length(pk);
%         pk = [pk; zeros(n,1)];
%         lk = [lk; zeros(n,1)];
%         
%         for jj = 0:n-1
%             IND = 117 - jj ;
%             [I,J] = ind2sub([3 39],IND);
%             isOrderComplete(I,J) = 0;
%         end
%         
%     end
        
    peak(:,ii)=pk;
    loc(:,ii)=lk;

end
% % 
% figure(1000);
% imagesc(flat)
% colormap gray;
% hold on;
% set(gca,'ydir','normal')

% fit the locations of each peak across an order 
for ii = 1:3*36
    ind = loc(ii,:) ~= 0;
   
    
%     plot(xx(ind),loc(ii,ind),'x','markersize',18)
    
    alpha=polyfit(xx(ind),loc(ii,ind),2);
    trace(:,ii) = polyval(alpha,1:4096);
%     plot(1:4096,trace(:,ii),'--','linewidth',2)

end

end