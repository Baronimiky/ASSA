function  u = input_f(t, amp,t_pulse)

if t>= t_pulse(1) && t<=t_pulse(2)
        force = amp;
 else 
     force = 0; 
 end
 u = @(t) force; 
end

