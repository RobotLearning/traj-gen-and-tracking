% basis functions are unscaled gaussians
function out = sqExp(x,h,c)
out = exp(-h * (x - c)^2);
end
