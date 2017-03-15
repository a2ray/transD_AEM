function P = SkyTEMDataProc(fileName,SoundingNum)

% This function processes the data from SkyTEM _dat.xyz files. The fileName
% specified in the input is the path to the _dat.xyz file. The SoundingNum
% is the line number corresponding to the particular sounding you want to
% invert. P is a structure that contains the data, the gate times, and the
% altitude of the loop for the specified sounding

%location of the _dat.xyz file to read data from
%fileName = 'TaylorGlacierSkyTEMdata/TaylorGlacier_dat.xyz';
%this strips away all the header information
!grep -v "/" TaylorGlacierSkyTEMdata/TaylorGlacier_dat.xyz > RawOut.txt
%this gets the gate times
!grep "/ " TaylorGlacierSkyTEMdata/TaylorGlacier_dat.xyz > tmptimes2.txt
!grep -v "LINE" tmptimes2.txt > tmptimes.txt

%load the data
load RawOut.txt
%reading in the gate times, with the leading '/' is trickier
fileID = fopen('tmptimes.txt');
tmp = textscan(fileID,'%s',1);  %first, read the useless leading '/'
tmp = textscan(fileID,'%f');    %now, read all the gate times following it
times = tmp{1};

%dBz/dt begins in column 18
RawData = RawOut(:,18:(18+37));
DataErr = RawOut(:,56:(56+37));
%Data type (HM or LM) indicator
DataType = RawOut(:,15);
FID = RawOut(:,4);

for j=1:size(RawData,1)
    for k=1:size(RawData,2)
        if( RawData(j,k) > 9000 )
            RawData(j,k) = NaN;
        end
        if( DataErr(j,k) > 9000 )
            DataErr(j,k) = NaN;
        end
    end
end

%let's plot up some data
% for j=1:11
%     loglog(times,RawData(j,:),'o')
%     hold on
% end
% xlabel('time (s)')
% ylabel('dBz/dt (T/s)')
% title('Some raw data from SkyTEM')

%remove temporary files
!rm tmptimes2.txt
!rm tmptimes.txt
!rm RawOut.txt

%write the extracted data and the gate times to file
%!rm GateTimes.txt
%!rm SkyTEMdata.txt
%dlmwrite('GateTimes.txt',times)
%dlmwrite('SkyTEMdata.txt',RawData)

if( DataType(SoundingNum) == 1 )
    P.dataLM = RawData(SoundingNum,:);
    P.dataErrLM = DataErr(SoundingNum,:)';
    if( DataType(SoundingNum+1) == 2 && FID(SoundingNum+1) == FID(SoundingNum) )
        P.dataHM = RawData(SoundingNum+1,:);
        P.dataErrHM = DataErr(SoundingNum+1,:)';
    elseif( DataType(SoundingNum-1) == 2 && FID(SoundingNum-1) == FID(SoundingNum) )
        P.dataHM = RawData(SoundingNum-1,:);
        P.dataErrHM = DataErr(SoundingNum-1,:)';
    else
        P.dataHM = [ ];
        P.dataErrHM = [ ];
        fprintf('There is no HM data for this sounding\n')
    end
elseif( DataType(SoundingNum) == 2 )
    P.dataHM = RawData(SoundingNum,:);
    P.dataErrHM = DataErr(SoundingNum,:)';
    if( DataType(SoundingNum+1) == 1 && FID(SoundingNum+1) == FID(SoundingNum) )
        P.dataLM = RawData(SoundingNum+1,:);
        P.dataErrLM = DataErr(SoundingNum+1,:)';
    elseif( DataType(SoundingNum-1) == 1 && FID(SoundingNum-1) == FID(SoundingNum) )
        P.dataLM = RawData(SoundingNum-1,:);
        P.dataErrLM = DataErr(SoundingNum-1,:)';
    else
        P.dataLM = [ ];
        P.dataErrLM = [ ];
        fprintf('There is no LM data for this sounding\n')
    end
else
    fprintf('Data error: data is neither HM or LM!\n')
end


% P.dataHM = RawData(SoundingNum,:)';   %High-mode data for this sounding
% P.dataLM = RawData(SoundingNum+1,:)';  %Low-mode data for this sounding
% P.dataErrHM = DataErr(SoundingNum,:)';
% P.dataErrLM = DataErr(SoundingNum+1,:)';
P.times = times;
P.alt = RawOut(SoundingNum,7);

P.XYcoord(:,1) = RawOut(:,2);
P.XYcoord(:,2) = RawOut(:,3);



