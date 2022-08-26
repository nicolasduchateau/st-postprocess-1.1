clc;
close all;

%% prepare data

ref = cell(numSHEET,1);
for xl=1:numSHEET
    xlREF = 1;
    sREF = 1;                      %%% one common reference, e.g. subject #1

    ref{xl}.index = sREF;
    ref{xl}.patient = patient{xlREF}{ref{xl}.index};
    ref{xl}.ECG_list = ref{xl}.patient.ECG_list;
    ref{xl}.numCP = ref{xl}.patient.STdata.numCP;
    
    f0 = ref{xl}.patient.manualECG(1)-ref{xl}.patient.manualECG(1)+1;   %%% frame beginning the cycle
    
    numSubjects = length(patient{xl});
    for s=[sREF,1:sREF-1,sREF+1:numSubjects]   %%% reference treated first
        if(patient{xl}{s}.toProcess == 1)
            disp([patient{xl}{s}.seriesIDlong,': time']);
            
%% Temporal interpolation
            patient{xl}{s}.STdata.numFramesCycle = patient{xl}{s}.manualECG(length(patient{xl}{s}.manualECG)) - patient{xl}{s}.manualECG(1) + 1;
            
            nTime_REF = ref{xl}.ECG_list(:,2);
            nTime_K = patient{xl}{s}.ECG_list(:,2);

            patient{xl}{s}.STdata.Txy_pos = NaN(2,size(ref{xl}.ECG_list,1),patient{xl}{s}.STdata.numCP);
            patient{xl}{s}.STdata.Tu_ro = NaN(2,2,size(ref{xl}.ECG_list,1),patient{xl}{s}.STdata.numCP);
            patient{xl}{s}.STdata.TPhi_ro = NaN(2,size(ref{xl}.ECG_list,1),patient{xl}{s}.STdata.numCP);
            patient{xl}{s}.STdata.Tv_ro = NaN(2,size(ref{xl}.ECG_list,1),patient{xl}{s}.STdata.numCP);
            patient{xl}{s}.STdata.TS_longRO = NaN(2,size(ref{xl}.ECG_list,1),patient{xl}{s}.STdata.numCP);

            for p=1:patient{xl}{s}.STdata.numCP
                for nDim=1:2
                    %% xy_pos
                    dataIN = squeeze( patient{xl}{s}.STdata.xy_pos(nDim,patient{xl}{s}.manualECG(1):patient{xl}{s}.manualECG(length(patient{xl}{s}.manualECG)),p) );
                    extrapVal = 0; removeLast = 0;
                    dataOUT = interp1([nTime_K(1:length(nTime_K)-removeLast)], dataIN(1:length(dataIN)-removeLast), [nTime_REF(1:length(nTime_REF)-removeLast)],'spline',extrapVal);
                    patient{xl}{s}.STdata.Txy_pos(nDim,:,p) = dataOUT;

                    %% u_ro
                    for nDim2=1:2
                        dataIN = squeeze( patient{xl}{s}.STdata.u_ro(nDim,nDim2,patient{xl}{s}.manualECG(1):patient{xl}{s}.manualECG(length(patient{xl}{s}.manualECG)),p) );
                        extrapVal = 1; removeLast = 0;
                        dataOUT = interp1([nTime_K(1:length(nTime_K)-removeLast)], dataIN(1:length(dataIN)-removeLast), [nTime_REF(1:length(nTime_REF)-removeLast)],'spline',extrapVal);
                        patient{xl}{s}.STdata.Tu_ro(nDim,nDim2,:,p) = dataOUT;
                    end
                    
                    %% Phi
                    dataIN = squeeze( patient{xl}{s}.STdata.Phi_ro(nDim,patient{xl}{s}.manualECG(1):patient{xl}{s}.manualECG(length(patient{xl}{s}.manualECG)),p) );
                    extrapVal = 0; removeLast = 0;
                    dataOUT = interp1([nTime_K(1:length(nTime_K)-removeLast)], dataIN(1:length(dataIN)-removeLast), [nTime_REF(1:length(nTime_REF)-removeLast)],'spline',extrapVal);
                    patient{xl}{s}.STdata.TPhi_ro(nDim,:,p) = dataOUT;                
                    
                    %% v
                    dataIN = squeeze( patient{xl}{s}.STdata.v_ro(nDim,patient{xl}{s}.manualECG(1):patient{xl}{s}.manualECG(length(patient{xl}{s}.manualECG)),p) )';
                    dataIN(length(dataIN)) = 0;
                    extrapVal = 0; removeLast = 0;
                    dataOUT = interp1([nTime_K(1:length(nTime_K)-removeLast)], dataIN(1:length(dataIN)-removeLast), [nTime_REF(1:length(nTime_REF)-removeLast)],'spline',extrapVal);
                    patient{xl}{s}.STdata.Tv_ro(nDim,:,p) = dataOUT;
                
                    %% S
                    dataIN = squeeze( patient{xl}{s}.STdata.S_longRO(nDim,patient{xl}{s}.manualECG(1):patient{xl}{s}.manualECG(length(patient{xl}{s}.manualECG)),p) );
                    extrapVal = 1; removeLast = 0;
                    dataOUT = interp1([nTime_K(1:length(nTime_K)-removeLast)], dataIN(1:length(dataIN)-removeLast), [nTime_REF(1:length(nTime_REF)-removeLast)],'spline',extrapVal);
                    patient{xl}{s}.STdata.TS_longRO(nDim,:,p) = dataOUT;
                end
            end
            ref{xl}.patient = patient{xlREF}{ref{xl}.index};

%% Spatial interpolation
            disp([patient{xl}{s}.seriesIDlong,': space']);
            
            CP_k = zeros(size(ref{xl}.ECG_list,1),patient{xl}{s}.STdata.numCP);
            CP_ref = zeros(size(ref{xl}.ECG_list,1),ref{xl}.numCP);
            for f=1:size(ref{xl}.ECG_list,1)
                tmp2 = squeeze( patient{xl}{s}.STdata.Txy_pos(:,f,1:patient{xl}{s}.STdata.numCP-1) );
                tmp1 = squeeze( patient{xl}{s}.STdata.Txy_pos(:,f,1:patient{xl}{s}.STdata.numCP) );
                tmp2 = [tmp2(:,1),tmp2];
                tmpD = (sum( (tmp1-tmp2).^2 , 1)).^(1/2);
                CP_k(f,:) = cumsum(tmpD);
                CP_k(f,:) = CP_k(f,:) ./ CP_k(f,size(CP_k,2));

                tmp2 = squeeze( ref{xl}.patient.STdata.Txy_pos(:,f0,1:ref{xl}.numCP-1) );
                tmp1 = squeeze( ref{xl}.patient.STdata.Txy_pos(:,f0,1:ref{xl}.numCP) );

                tmp2 = [tmp2(:,1),tmp2];
                tmpD = (sum( (tmp1-tmp2).^2 , 1)).^(1/2);

                CP_ref(f,:) = cumsum(tmpD);
                CP_ref(f,:) = CP_ref(f,:) ./ CP_ref(f,size(CP_ref,2));
            end

            patient{xl}{s}.STdata.RSTxy_pos = NaN(2,size(ref{xl}.ECG_list,1),ref{xl}.numCP);
            patient{xl}{s}.STdata.RSTu_ro = NaN(2,2,size(ref{xl}.ECG_list,1),ref{xl}.numCP);
            patient{xl}{s}.STdata.RSTPhi_ro = NaN(2,size(ref{xl}.ECG_list,1),ref{xl}.numCP);
            patient{xl}{s}.STdata.RSTv_ro = NaN(2,size(ref{xl}.ECG_list,1),ref{xl}.numCP);
            patient{xl}{s}.STdata.RSTS_longRO = NaN(2,size(ref{xl}.ECG_list,1),ref{xl}.numCP);
            
            for f=1:size(ref{xl}.ECG_list,1)
                for nDim=1:2
                    %% xy_pos
                    dataIN = squeeze( patient{xl}{s}.STdata.Txy_pos(nDim,f,:) );
                    extrapVal = 0; removeLast = 0;
                    dataOUT = interp1(CP_k(f,:),dataIN,CP_ref(f,:),'spline',extrapVal);
                    patient{xl}{s}.STdata.RSTxy_pos(nDim,f,:) = dataOUT;

                    %% u_ro
                    for nDim2=1:2
                        dataIN = squeeze( patient{xl}{s}.STdata.Tu_ro(nDim,nDim2,f,:) );
                        extrapVal = 1; removeLast = 0;
                        dataOUT = interp1(CP_k(f,:),dataIN,CP_ref(f,:),'spline',extrapVal);
                        patient{xl}{s}.STdata.RSTu_ro(nDim,nDim2,f,:) = dataOUT;
                    end
                    
                    %% Phi
                    dataIN = squeeze( patient{xl}{s}.STdata.TPhi_ro(nDim,f,:) );
                    extrapVal = 0; removeLast = 0;
                    dataOUT = interp1(CP_k(f,:),dataIN,CP_ref(f,:),'spline',extrapVal);
                    patient{xl}{s}.STdata.RSTPhi_ro(nDim,f,:) = dataOUT;
                    
                    %% v
                    dataIN = squeeze( patient{xl}{s}.STdata.Tv_ro(nDim,f,:) );
                    extrapVal = 0; removeLast = 0;
                    dataOUT = interp1(CP_k(f,:),dataIN,CP_ref(f,:),'spline',extrapVal);
                    patient{xl}{s}.STdata.RSTv_ro(nDim,f,:) = dataOUT;
                
                    %% S
                    dataIN = squeeze( patient{xl}{s}.STdata.TS_longRO(nDim,f,:) );
                    extrapVal = 1; removeLast = 0;
                    dataOUT = interp1(CP_k(f,:),dataIN,CP_ref(f,:),'spline',extrapVal);
                    patient{xl}{s}.STdata.RSTS_longRO(nDim,f,:) = dataOUT;
                end
            end
            clear patient{xl}{s}.STdata.Txy_pos;
            clear patient{xl}{s}.STdata.Tu_ro;
            clear patient{xl}{s}.STdata.TPhi_ro;
            clear patient{xl}{s}.STdata.Tv_ro;
            clear patient{xl}{s}.STdata.TS_longRO;
        end
    end
end   
    
