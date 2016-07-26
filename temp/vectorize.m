function val = vectorize(func,x)

fmat = func(x);
val = fmat(:);
end