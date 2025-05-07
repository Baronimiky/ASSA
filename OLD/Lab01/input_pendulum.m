function  u = input_pendulum(t, amp, t0)
 
if t>= t0
     force = amp;
 else 
     force = 0; 
 end
 u = @(t) force; 
end

