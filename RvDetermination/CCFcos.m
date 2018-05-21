function [ccfresult]=CCFcos(m_l, m_h, wav, spec, v_r,weight)

%     this function calculates the cross-correlation function;
%     of a spectrum(spe, at wavlenths wav) with a mask. m_l;
%     and m_h contain the beginning and eml of regions where the mask;
%     equals 1, ordered.;
% snw =1;
sn = ones(size(wav));
snw = 1;
% weight = ones(size(m_l));
ccfresult=0;
%     speed of light, km./s;
c = 2.99792458e5;

%     doppler factor;
gamma = sqrt(1. + v_r ./ c) ./ sqrt(1. - v_r ./ c);

%     doppler shift mask;

m_lloc = m_l.* gamma;
m_hloc = m_h.* gamma;

% hold on
% for ii = 1:1
% plot(m_lloc(ii).*[1 1], [0 3])
% plot(m_hloc(ii).*[1 1], [0 3])
% end


% HW = 0.5*(m_hloc - m_lloc);

%     i marks where we are in terms of masks;
ii = 1;

cond = 0;
N = length(m_lloc);

% pix_init(jj)= zeros(length(wav)-1,1);
% pix_end(jj)= zeros(length(wav)-1,1);



for jj = 2:length(wav)-1


%     fprintf('pix_init: %.4f pix_end: %.4f\n',pix_init(jj),pix_end(jj))
    
end    

for jj = 2:length(wav)-1
    pix_init = 0.5.*(wav(jj-1) + wav(jj));
    pix_end  = 0.5.*(wav(jj) + wav(jj+1));
        
    while (m_hloc(ii) < pix_init) && (cond == 0)
        if(ii == N)
            cond = 1;
        end
        
        if(cond == 0)
            ii = ii + 1;
        end
    end


%     NORM = HW(I) 
%     ARG  = PI / (2.0D0*HW(I))

    if((pix_end < m_hloc(ii)) && (pix_init > m_lloc(ii))) 

        ccfresult = ccfresult + spec(jj) .* weight(ii) .* sn(jj);
        snw = snw + sn(jj).*weight(ii);
    elseif((pix_end < m_hloc(ii)) && (pix_init< m_lloc(ii)) && (pix_end> m_lloc(ii))) 

            fraction_ml =(pix_end- m_lloc(ii)) ./(pix_end- pix_init);
            ccfresult = ccfresult + spec(jj) .* weight(ii) .* fraction_ml .* sn(jj);
            snw = snw + fraction_ml.*sn(jj).*weight(ii);
    elseif((pix_end > m_hloc(ii)) && (pix_init > m_lloc(ii)) && (pix_init < m_hloc(ii))) 

            fraction_ml =(m_hloc(ii) - pix_init) ./(pix_end - pix_init);
            ccfresult = ccfresult + spec(jj) .* weight(ii) .* fraction_ml .* sn(jj);
            snw = snw + fraction_ml.*sn(jj).*weight(ii);
    elseif((pix_end > m_hloc(ii)) && (pix_init < m_lloc(ii)))

             fraction_ml =(m_hloc(ii) - m_lloc(ii)) ./(pix_end- pix_init);
             ccfresult = ccfresult + spec(jj) .* weight(ii) .* fraction_ml .* sn(jj);
             snw = snw + fraction_ml.*sn(jj).*weight(ii);
    end 
end
%      ccfresult = ccfresult ./ snw;

end