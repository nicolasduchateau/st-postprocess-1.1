function         ECG_list = getECGnormalizedScale(manualECG)

    manualECG = manualECG';
    startingImageNumber = manualECG(1);
    endingImageNumber = manualECG(length(manualECG));       
    numFrames = endingImageNumber - startingImageNumber + 1;

    ECG_list = zeros(numFrames,2);
    ECG_list(:,1) = (1:numFrames) + startingImageNumber - 1;

    ECGKeyPoints = manualECG;

    k = 1;
    for f = 1:numFrames
        if ( ( ECG_list(f,1) > ECGKeyPoints(k+1,1) ) && ( ECGKeyPoints(k+1,1) < ECGKeyPoints(length(ECGKeyPoints),1) ) )
            k = k + 1;
        end
        ECG_list(f,2) = (ECG_list(f,1) - ECGKeyPoints(k,1)) / ...
                        (ECGKeyPoints(k+1,1) - ECGKeyPoints(k,1)) + ...
                        (k-1);
    end

end