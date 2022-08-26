close all;
clc;

CPtoPlot = 10;

for xl=1:numSHEET
    framesToPlot = ( ref{xl}.patient.manualECG(1):ref{xl}.patient.manualECG(6) ) - ref{xl}.patient.manualECG(1) + 1;
    
    %% plot data at a given location
    subplot(1,numSHEET,xl);
    numSubjects = length(patient{xl});
    for s=1:numSubjects
        if(patient{xl}{s}.toProcess == 1)
            nDim = 2;   %%% plot 2nd component = longitudinal
            dataToPlot = patient{xl}{s}.STdata.RSTPhi_ro;   %%% displacement in radial/longitudinal coordinates
            dataToPlot = dataToPlot(nDim,framesToPlot,CPtoPlot);
            plot( dataToPlot ); hold on;
        end
    end
    
    drawnow;
    axis manual; axis square; 
    for e=1:length(ref{xl}.patient.manualECG);
        idx = ref{xl}.patient.manualECG(e) - ref{xl}.patient.manualECG(1) + 1;
        plot([idx,idx],[-1000,1000],'k');
    end
    axis01 = axis;
    axis([1,length(framesToPlot),axis01(3),axis01(4)]);

end
