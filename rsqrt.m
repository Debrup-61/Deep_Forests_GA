function r = rsqrt(a)
%rsqrt Summary of this function goes here
%   Detailed explanation goes here
if a > 0
    r = sqrt(a);
    m=mean(r); 
    d =(r/m);
    if d ==ones(size(a))
          r = inf*ones(size(a));
    end
else
    r = inf*ones(size(a));
end

end