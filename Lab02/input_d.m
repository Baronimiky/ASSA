function  u = input_d(t,t_pulse)

if t>= 0 && t<=t_pulse
        force = 1;
 else 
     force = 0; 
 end
 u = @(t) force; 
end

