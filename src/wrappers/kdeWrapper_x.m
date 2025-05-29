function u = kdeWrapper_x(a,kernelType,probType,bandwidthType,SupportType,nPoints)
      [~,u] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType,NumPoints=nPoints);
end