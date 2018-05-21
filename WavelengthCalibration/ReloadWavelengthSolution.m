function [wave_cal] = ReloadWavelengthSolution(coeffs,trace)

load(coeffs); 
cfs = wave_coeff(:,:,trace); %trace num
xxf=((0.5:4095.5)-2048)/100;

for ii = 1:36
    zem(:,ii) = 1e4*polyval(cfs(ii,:),xxf);
end

wave_cal = zem;

end