function r = sine(a)
%SIN Summary of this function goes here
%   Detailed explanation goes here
   r=sin(a);
   if r == zeros(size(a))
     r = inf*ones(size(a));
   else
     m=mean(r); 
     d =(r/m);
     
     if d ==ones(size(a))
          r = inf*ones(size(a));
     end
   
   
   end
end   
   

