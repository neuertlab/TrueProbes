function u = kdeWrapper_y(a,kernelType,probType,bandwidthType,SupportType)
      [u,~] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType);
end