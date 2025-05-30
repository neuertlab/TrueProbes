function [paramsets] = sample_parameters(muu,covv,nSamples,constrains,posterior,lin_dist_slope,sampling,resample_dist,re_sampling_outliers)
    UB=constrains.log10UB; LB=constrains.log10LB; RANG=UB-LB; MED=median([UB; LB],1); delta=constrains.delta*RANG; 
        % sample parameters from multivariate normal random  number distribution
    switch sampling
    case "mvnrnd" 
        paramsets = mvnrnd(muu, covv, nSamples);
    case "disperse" 
        paramsets = disperse(resample_dist,posterior,delta);
    end
%     [fig] = plot_MVN(paramsets,[1:4],constrains,1,dirName);
%     FIGsFrames(1)=getframe(fig); 
    for j=1
    paramsets0=paramsets; 
    for i=1:length(muu)
        this_parameter = paramsets0(:,i); %this_parameter
        % parameters at bounds proximity (within 25% from the bounds) or outside                
            % set them to randn
    %     INDEX_U=find(this_parameter>(UB-(1/8)*RANG)); 
    %     INDEX_L=find(this_parameter<(LB+(1/8)*RANG)); 
    %     this_parameter(INDEX_U)=(1/3)*.5*RANG(1)*randn(length(INDEX_U),1);
    %     this_parameter(INDEX_L)=(1/3)*.5*RANG(1)*randn(length(INDEX_L),1); 
    %     [fig] = plot_MVN(this_parameter,muu,covv,[1:4],constrains,4,dirName);

            % set them to randn probabilistically (closer to the bounds, more likely to be resampled)           
        INDEX_U=find(this_parameter>(UB(i)-(1/4)*RANG(i))); % ID for parameters too big p>c
        INDEX_L=find(this_parameter<(LB(i)+(1/4)*RANG(i))); % ID for parameters too small p<c

        xu=abs((1/3)*(1/8)*RANG(i)*randn(length(INDEX_U),1)); % normal random range from UB
        xl=abs((1/3)*(1/8)*RANG(i)*randn(length(INDEX_L),1)); % normal random range from LB
        
%         xu=abs((1/3)*(1/4)*RANG(i)*randn(length(INDEX_U),1)); % normal random range from UB
%         xl=abs((1/3)*(1/4)*RANG(i)*randn(length(INDEX_L),1)); % normal random range from LB

        pp_u=this_parameter(INDEX_U); % parameters too big
        pp_l=this_parameter(INDEX_L); % parameters too small

        INDEX_UU = find(pp_u>(UB(i)-xu)); % ID for parameters too big p>x
        INDEX_LL = find(pp_l<(LB(i)+xl)); % ID for parameters too small p<x
        
        switch sampling            
        case "mvnrnd" % resample lognormal around 0
            if (re_sampling_outliers=="lognormal_at_mu" & LB(i)<muu(i) & muu(i)<UB(i) )
                pp_u(INDEX_UU)=muu(i) + 1e-5*randn(length(INDEX_UU),1); % resample at mu
                pp_l(INDEX_LL)=muu(i) + 1e-5*randn(length(INDEX_LL),1); % resample at mu
            else
                pp_u(INDEX_UU)=(1/4)*.5*RANG(i)*randn(length(INDEX_UU),1); % resample at MED
                pp_l(INDEX_LL)=(1/4)*.5*RANG(i)*randn(length(INDEX_LL),1); % resample at MED            
            end
        case "disperse" 
            switch re_sampling_outliers
            case "uniform_each_half"
            % re-sample uniform at - or + half range depending on which side this parameter sticks out
            pp_u(INDEX_UU)= MED(i) + .49*RANG(i)*randl(0,[length(INDEX_UU) 1]); % resample  (.49 instead of .5) 2% off the bounds
            pp_l(INDEX_LL)= MED(i) - .49*RANG(i) + .49*RANG(i)*randl(0,[length(INDEX_LL) 1]); % resample
            case "linear_each_half"
            % re-sample linear at - or + half range increasing toward L or U depending on which side this parameter sticks out
            pp_u(INDEX_UU)= MED(i) + .49*RANG(i)*randl(lin_dist_slope,[length(INDEX_UU) 1]); % resample  (.49 instead of .5) 2% off the bounds
            pp_l(INDEX_LL)= MED(i) - .49*RANG(i) + .49*RANG(i)*randl(-lin_dist_slope,[length(INDEX_LL) 1]); % resample
            case  "lognormal_at_mu"
            pp_u(INDEX_UU)= posterior.PARAS(1,i) + 1e-5*randn(length(INDEX_UU),1); % resample at the best pars (i-1)
            pp_l(INDEX_LL)= posterior.PARAS(1,i) + 1e-5*randn(length(INDEX_LL),1); % resample at the best pars (i-1)
            end
        end

        this_parameter(INDEX_U) = pp_u; % replace
        this_parameter(INDEX_L) = pp_l; % replace

        paramsets(:,i) = this_parameter;
    end
%     [fig] = plot_MVN(paramsets,[1:4],constrains,j,dirName);
%     FIGsFrames(j)=getframe(fig); 
    end
%     implay(FIGsFrames, 1);
%     VideoName = [dirName, 'parameters_15iteration_resampling'];
%     v = VideoWriter(VideoName,'MPEG-4'); 
%     v.FrameRate = 0.5;
%     open(v)   
%     writeVideo(v,FIGsFrames);
%     close(v); 
end
