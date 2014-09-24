% basis functions are unscaled gaussians
function out = vonMises(phi,h,c)
out = exp(h * (cos(phi - c) - 1));
end