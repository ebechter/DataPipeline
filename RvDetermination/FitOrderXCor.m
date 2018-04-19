function [order_RV] = FitOrderXCor(start_order,spec_order,XCvelocities,vel0,XCF,masksize,plotflag,injected)



for kk = start_order:spec_order
    
    gstartpoints = [200 vel0 1 0];
    [~,order_RV(kk)] = GaussianFit(XCvelocities,XCF(kk,:)',gstartpoints,plotflag,'k',1);
    
    fprintf('%i Individual RV (m/s): %.5f mask size: %i \n',kk, order_RV(kk)*1000-injected,masksize(kk))
    
end

    fprintf('Mean (m/s): %.5f\n',mean(order_RV)*1000)
