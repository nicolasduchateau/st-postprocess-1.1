function STdata_OUT = getSTdataXY(STdata_IN,manualECG,options)

numFrames = STdata_IN.numFrames;
numCP = STdata_IN.numCP;
xy_pos = STdata_IN.xy_pos;              %% position in Cartesian system of coordinates

%% output data
u_ro = NaN(2,2,numFrames,numCP);        %% local system of coordinates
xy_posOUT = NaN(2,numFrames,numCP);     %% position in Cartesian system of coordinates
Phi_ro = NaN(2,numFrames,numCP);        %% displacement in Cartesian system of coordinates
v_ro = NaN(2,numFrames,numCP);          %% velocity in Cartesian system of coordinates
S_longRO = NaN(2,numFrames,numCP);      %% strain along the myocardium in Cartesian system of coordinates

timeInterval = ( STdata_IN.endTime - STdata_IN.beginTime ) / ( manualECG(length(manualECG)) - manualECG(1) );

%% Drift correction
if options.correctDrift == 1
    xy_pos_old = xy_pos;
    for f=1:numFrames
        for p=1:numCP
            xy_pos(:,p,f) = xy_pos(:,p,f) - ( xy_pos_old(:,p,manualECG(length(manualECG))) - xy_pos_old(:,p,manualECG(1)) ) * (f - manualECG(1)) / (manualECG(length(manualECG)) - manualECG(1));
        end
    end
end

%% Correction according to center of mass
if options.correctCenter == 1
    xy_pos_old = xy_pos;
    for f=1:numFrames
        if STdata_IN.is4CH == 1         %% for SAX data only
            numCP2 = 2*floor(numCP/2);
            if numCP == numCP2
                idx1 = 1:floor(numCP/2);
                idx2 = numCP:-1:floor(numCP/2)+1;
            else
                idx1 = [1:floor(numCP/2),floor(numCP/2)+1];
                idx2 = [numCP:-1:floor(numCP/2)+1];
            end            
            centerLine = ( xy_pos(:,idx1,f) + xy_pos(:,idx2,f) )/2;
            corrT = centerLine;
%             corrT(2,:) = centerLine(2,length(centerLine));
            corrT(2,:) = 0;
            xy_pos(:,idx1,f) = xy_pos_old(:,idx1,f) - corrT;
            xy_pos(:,idx2,f) = xy_pos_old(:,idx2,f) - corrT;
        else
            barycenter_f = mean( xy_pos(:,:,f) , 2 );
            xy_pos(:,:,f) = xy_pos_old(:,:,f) - repmat(barycenter_f,1,numCP);
        end
    end
end

%% Get radial / orthogonal (longitudinal or circumferential) data
for f=1:numFrames
    for p=1:numCP
        if STdata_IN.is4CH == 0         %% for SAX data only
            switch p
                case 1                  %% first data point
                    uo = (xy_pos(:,p+1,f) - xy_pos(:,numCP,f))/2;
                case numCP              %% last data point
                    uo = (xy_pos(:,1,f) - xy_pos(:,p-1,f))/2;
                otherwise
                    uo = (xy_pos(:,p+1,f) - xy_pos(:,p-1,f))/2;
            end
        else                            %% for 4CH data only
            switch p
                case 1                  %% first data point
                    uo = xy_pos(:,p+1,f) - xy_pos(:,p,f);
                case numCP              %% last data point
                    uo = xy_pos(:,p,f) - xy_pos(:,p-1,f);
                otherwise
                    uo = (xy_pos(:,p+1,f) - xy_pos(:,p-1,f))/2;
            end
        end
        uo = uo ./ norm(uo);
        Rm = [0,-1;1,0];
        ur = Rm * uo;
        if ((STdata_IN.Flip == 0)&&(STdata_IN.is4CH == 0))
            ur = -ur;
        end
        if STdata_IN.is4CH == 1
            if p>floor(numCP/2)
                uo = -uo;
            end
        end
        u_ro(:,:,f,p) = [ur,uo];
    end
end

%% displacement in Cartesian system of coordinates
for f=1:numFrames
    for p=1:numCP
        uxy = xy_pos(:,p,f) - xy_pos(:,p,manualECG(1));
        Phi_ro(:,f,p) = inv(u_ro(:,:,f,p)) * uxy;
    end
end

%% velocity in Cartesian system of coordinates
for f=1:numFrames-1
    for p=1:numCP
        vxy = ( xy_pos(:,p,f+1) - xy_pos(:,p,f) ) / timeInterval;
        v_ro(:,f,p) = inv(u_ro(:,:,f,p)) * vxy;
    end
end

%% strain in Cartesian system of coordinates
for f=1:numFrames
    for p=1:numCP
        if STdata_IN.is4CH == 0         %% for SAX data only
            switch p
                case 1
                    tmp_f = (xy_pos(:,p+1,f) - xy_pos(:,numCP,f))/2;
                    tmp_1 = (xy_pos(:,p+1,manualECG(1)) - xy_pos(:,numCP,manualECG(1)))/2;
                case numCP
                    tmp_f = (xy_pos(:,1,f) - xy_pos(:,p-1,f))/2;
                    tmp_1 = (xy_pos(:,1,manualECG(1)) - xy_pos(:,p-1,manualECG(1)))/2;
                otherwise
                    tmp_f = ( xy_pos(:,p+1,f) - xy_pos(:,p-1,f) )/2;
                    tmp_1 = ( xy_pos(:,p+1,manualECG(1)) - xy_pos(:,p-1,manualECG(1)) )/2;
            end
        else                            %% for 4CH data only
            switch p
                case 1
                    tmp_f = (xy_pos(:,p+1,f) - xy_pos(:,p,f))/2;
                    tmp_1 = (xy_pos(:,p+1,manualECG(1)) - xy_pos(:,p,manualECG(1)))/2;
                case numCP
                    tmp_f = (xy_pos(:,p,f) - xy_pos(:,p-1,f))/2;
                    tmp_1 = (xy_pos(:,p,manualECG(1)) - xy_pos(:,p-1,manualECG(1)))/2;
                otherwise
                    tmp_f = ( xy_pos(:,p+1,f) - xy_pos(:,p-1,f) )/2;
                    tmp_1 = ( xy_pos(:,p+1,manualECG(1)) - xy_pos(:,p-1,manualECG(1)) )/2;
            end
        end
        %% from Cartesian to local u_ro system of coordinates
        tmpINV = inv(u_ro(:,:,f,p));
        tmp_f = tmpINV * tmp_f;
        tmp_1 = tmpINV * tmp_1;
        tmpS = tmp_f ./ tmp_1;
        tmpS(1) = 0;            %% 1D data only
        %% from local u_ro to Cartesian system of coordinates
        S_longRO(:,f,p) = tmpS;
    end
end

%% Re-indexing position in Cartesian system of coordinates
for f=1:numFrames
    for p=1:numCP
        xy_posOUT(:,f,p) = xy_pos(:,p,f);
    end
end

%% Storing data
STdata_OUT = STdata_IN;
STdata_OUT.timeInterval = timeInterval;
STdata_OUT.xy_pos = xy_posOUT;
STdata_OUT.u_ro = u_ro;
STdata_OUT.Phi_ro = Phi_ro;
STdata_OUT.v_ro = v_ro;
STdata_OUT.S_longRO = S_longRO;

end