function [sel, ptheta] = lojResultant(resp, f, thvals_deg)
%% Get OSI

    % resp: response vector

    % resultant
    result  = sum(resp.*exp(1i*f*thvals_deg*pi/180)); 
    
    % selectivity index
    sel     = abs(result)/sum(abs(resp)); 
    
    % est preffered ori
    ptheta    = rad2deg(angle(result)/f); 
    ptheta(ptheta < 0) = ptheta(ptheta < 0) +  360/f;

end