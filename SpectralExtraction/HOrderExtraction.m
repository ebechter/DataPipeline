function [sumorder] = HOrderExtraction(science,window,norm)

%% Sum over columns technique 

order = zeros(window*2+1,4096);

[pk,lk] = findpeaks(sum(science,2));

start = lk;





% loop through all the orders
for jj = 1:36
    
    for ii = 1:4096
        
%         order(:,ii)=science(round(trace(ii,jj))-window:round(trace(ii,jj))+window,ii); % collect flux in extraction window
        order(:,ii)=science(round(start(jj))-window:round(start(jj))+window,ii); % collect flux in extraction window

    end
    
    if norm ==1
        sumorder(:,jj)=sum(order)./max(sum(order)); % combine window
    else
        sumorder(:,jj)=sum(order);
    end
end