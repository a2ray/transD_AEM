function P = SkyTEMInvProc(fileName,FIDnum)

%this strips away all the header information
str = ['grep -v "/" ' fileName '_inv.xyz > RawOut.txt'];
system(str);

%load everything in the file
load RawOut.txt

%the FID numbers of the soundings
FID = RawOut(:,4);

%pull out the smoothed inversion result for this FIDnum
foundit = false;
for l=1:length(FID)
    if( FID(l) == FIDnum )
        fprintf('found FID number %d %d %d \n\n',l,FID(l),FIDnum)
        P.inv = RawOut(l,21:49);
        P.z = RawOut(l,108:136);
        P.zbot = RawOut(l,137:164);
        P.thick = RawOut(l,165:192);
        P.DOI = RawOut(l,end-1:end);
        foundit = true;
    end
end

if( foundit == false )
    fprintf('failed to find FID number %d \n\n',FIDnum)
end
