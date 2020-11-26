function [atws] = ScaleMasconToGrid(mgtws,mscal)
% function for simply scaling GRACE TWS into EALCO Model grid
%
% Input:
%   mgtws -- the GRACE TWS value
%   mscal -- scale factor for leakage error
%   
% Output:
%   atws  -- the scaled TWS data 
%

atws = mscal*mgtws;
    

