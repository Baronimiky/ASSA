function  u = input_pendulum_impulse(t, amp, t0)
    % 
    % if  t==t0
    %      force = amp;
    %  else 
    %      force = 0; 
    %  end
    %  u = force; 

    u = @(t) (t < t0 ) *amp; % Impulso breve di durata ~1ms
     
end

