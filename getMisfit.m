
%Large Loop misfit function
function misfit = getMisfit(x,S)

misfit = zeros(1,2);%for chi2 and RMS

%data usually ordered as [Er,Hb,Ez,Eb,Hr,Hz];
data = S.data;

% Get Bz for loop source:
Bz = get_field(S,x);
 
error_bar = S.sd;
 
m = (S.data - Bz)./S.sd;
 

misfit = [m'*m, sqrt(m'*m/length(S.data))];

end