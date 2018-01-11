function [sumorder] = BasicOrderExtraction(science,trace,window,norm)

%% Sum over columns technique 
order = zeros(window*2+1,4096);

% loop through all the orders
for jj = 1:36
   
    for ii = 1:4096
        
        order(:,ii)=science(round(trace(ii,jj))-window:round(trace(ii,jj))+window,ii); % collect flux in extraction window
        
    end
    
    if norm ==1
        sumorder(:,jj)=sum(order)./max(sum(order)); % combine window
    else
        sumorder(:,jj)=sum(order);
    end
end