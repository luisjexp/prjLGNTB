function [P, Kbin] = KernelRegprops(K)

    nrm = @(x) ( x-min(x(:)) )/range(x(:));
    Kbin = nrm(K{:})>.85;
    P = regionprops(nrm(K{:})>.85, 'Area', 'centroid', 'BoundingBox', 'Image', 'PixelIdxList', 'PixelList');
    
    P = P( [P.Area] == max([P.Area]) );

    if isempty(P)
        P = struct('Area', nan, 'Centroid', [nan nan], 'BoundingBox', nan(1,4),  'Image', [], 'PixelIdxList', [], 'PixelList', []);        
    end
    
end