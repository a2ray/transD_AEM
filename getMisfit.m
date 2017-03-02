
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
    
%     %temporary plotting code
%     eHM = S.HighMode.data+S.HighMode.sd;
%     eHM = [ eHM ; S.HighMode.data-S.HighMode.sd];
%     eLM = S.LowMode.data+S.LowMode.sd;
%     eLM = [ eLM ; S.LowMode.data-S.LowMode.sd];
%     clf
%     figure
%     loglog(S.HighMode.times,S.HighMode.data,'ro')
%     hold on
%     loglog(S.LowMode.times,S.LowMode.data,'bo')
%     loglog([S.HighMode.times;S.HighMode.times],eHM,'-')
%     loglog([S.LowMode.times;S.LowMode.times],eLM,'-')
%     loglog(S.HighMode.times,abs(Bz(1:19)),'*')
%     loglog(S.LowMode.times,abs(Bz(20:end)),'*')
%     hold off
   
end

%     m = (abs(data(:)) - abs(Bz(:)))./(error_bar(:));
%     misfit = [m'*m, sqrt(m'*m/length(error_bar))];
    residHM = (abs(S.HighMode.data)-abs(Bz(1:19)'))./S.HighMode.sd;
    residLM = (abs(S.LowMode.data)-abs(Bz(20:end)'))./S.LowMode.sd;
    rms_residHM = sqrt(residHM*residHM'/length(residHM));
    rms_residLM = sqrt(residLM*residLM'/length(residLM));
    rms_total = sqrt((residHM*residHM' + residLM*residLM')/length(residHM) + length(residLM));
    misfit = [ (rms_residHM*rms_residHM' + rms_residLM*rms_residLM')/2 , rms_total ];
    %fprintf('rmsHM, rmsLM, rmsTot: %f, %f, %f\n',rms_residHM,rms_residLM,rms_total)
    %keyboard
    
end



