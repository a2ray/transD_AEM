
%Large Loop misfit function
function misfit = getMisfit(x,S)
 
% Get Bz for loop source:
Bz = get_field(S,x);
 
if isfield(S,'data')  % Simple format
    
    data = S.data;
    error_bar = 2*S.sd;

else % SkyTEM High and Low Mode data:
    
    data = [S.HighMode.data S.LowMode.data];
    error_bar = 2*[S.HighMode.sd S.LowMode.sd]; 
   
end

    m = (abs(data(:)) - abs(Bz(:)))./(error_bar(:));
    misfit = [m'*m, sqrt(m'*m/length(error_bar))];
    
end