function milliPause(time)
%   MILLIPAUSE(X) pauses for X milliseconds.
%   
%   The original author of this should be attributed to  Marcel Kraemer.
%
%   Example of use:
%       %Pause for 100 milliseconds
%       milliPause(100); 
%
    t = time/1000;
    tic;
    while toc < t
    end
end
