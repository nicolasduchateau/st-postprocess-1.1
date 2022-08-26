clc;
clear all;
close all;

%% prepare data
parentFolder = cd(cd('..'));
folderRoot = [parentFolder,'/Data/'];   %%% Folder where the datafile is located
ExcelFileName = [folderRoot,'DEMO_DATA.xls'];   %%% datafile name

[typ,sheetsTitle] = xlsfinfo(ExcelFileName);
numSHEET = length(sheetsTitle);   %%% each sheet may correspond to a different protocol or stage (e.g. BASAL, LBBB, CRT, or VOL / PIG)
patient = cell(numSHEET,1);
options.lengthECG = 6;          %%% number of ECG events = 6 (Q1, MVC, AVO, AVC, MVO, Q2)
options.correctDrift = 1;       %%% use of drift correction (recommended)
options.correctCenter = 0;      %%% use of position correction wrt center of mass (not recommended for 4CH data)

for xl=1:numSHEET
    [num,txt,raw] = xlsread(ExcelFileName,sheetsTitle{xl});
    numSubjects = size(raw,1)-1;
    patient{xl} = cell(numSubjects,1);
    for s=1:numSubjects
%% reading Excel file information column per column
        i=1;                    
        patient{xl}{s}.population = raw{s+1,i}; i=i+1;
        patient{xl}{s}.seriesNumber = str2double(raw{s+1,i}); i=i+1;
        patient{xl}{s}.status = raw{s+1,i}; i=i+1;
        patient{xl}{s}.view = raw{s+1,i}; i=i+1;
        for e=1:options.lengthECG
            patient{xl}{s}.manualECG(e) = str2double(raw{s+1,i}) + 1; i=i+1;
        end
        STdata.is4CH = str2double(raw{s+1,i}); i=i+1;
        STdata.Flip = str2double(raw{s+1,i}); i=i+1;
        
        %%% (you can check the folder organization in the example provided)
        patient{xl}{s}.seriesIDshort = [patient{xl}{s}.population,get_extension(patient{xl}{s}.seriesNumber)];
        patient{xl}{s}.seriesID = [patient{xl}{s}.seriesIDshort ,'_',patient{xl}{s}.status];
        patient{xl}{s}.seriesIDlong = [patient{xl}{s}.seriesID,'_',patient{xl}{s}.view];
        patient{xl}{s}.pathE = [folderRoot,patient{xl}{s}.seriesIDshort,'/'];
    
%% reading speckle tracking exported data (CSV file)
        filename_read = [patient{xl}{s}.pathE,patient{xl}{s}.seriesIDlong,'.CSV'];

        if(exist(filename_read,'file'))
            disp(patient{xl}{s}.seriesIDlong);
            patient{xl}{s}.toProcess = 1;
            fid = fopen(filename_read);
            scan_format = '%s';
            dataTMP00 = textscan(fid,'%s',7,'Delimiter','\r','emptyvalue',NaN);
            dataTMP01 = textscan(fid,'%[^=]%1s%n%[^=]%1s%n%[^=]%1s%n%[^=]%1s%n',1,'Delimiter','\r','emptyvalue',NaN);
            dataTMP02 = textscan(fid,'%s',2,'Delimiter','\r','emptyvalue',NaN);
            dataTMP1 = textscan(fid,'%f',2,'Delimiter','\t','emptyvalue',NaN);
            dataTMP02 = textscan(fid,'%s',1,'Delimiter','\r','emptyvalue',NaN);
            dataTMP2 = textscan(fid,'%f','Delimiter',',','emptyvalue',NaN);

            fclose(fid);
            
            STdata.FR = dataTMP01{3};
            STdata.beginTime = dataTMP01{6};
            STdata.endTime = dataTMP01{9};
            STdata.ESTime = dataTMP01{12};
            
            STdata.numFrames = dataTMP1{1}(1);
            STdata.numCP = dataTMP1{1}(2);
            STdata.xy_pos = reshape(dataTMP2{1},2,STdata.numCP,STdata.numFrames);
            if (STdata.is4CH == 0)        %%% for SAX data only
                STdata.numCP = STdata.numCP - 1;
                STdata.xy_pos = STdata.xy_pos(:,1:STdata.numCP,:);
            end
    
%% get Phi / v / S / SR data
            STdata = getSTdataXY(STdata,patient{xl}{s}.manualECG,options);   %%% calling the function that will process the data
            patient{xl}{s}.STdata = STdata;

            patient{xl}{s}.ECG_list = getECGnormalizedScale(patient{xl}{s}.manualECG);   %%% normalization of ECG scale (not yet the data)
        else
            patient{xl}{s}.toProcess = 0;
        end
    end
end

