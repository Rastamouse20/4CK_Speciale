function [xfilt] = filterApply(x,filter)

    xfilt = (filter.H*x)./filter.Hs;
    
end