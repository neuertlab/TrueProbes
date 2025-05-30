<<<<<<< HEAD
<<<<<<< HEAD
function u = kdeWrapper_x(a,kernelType,probType,bandwidthType,SupportType,nPoints)
      [~,u] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType,NumPoints=nPoints);
=======
function u = kdeWrapper_x(a,kernelType,probType,bandwidthType,SupportType,nPoints)
      [~,u] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType,NumPoints=nPoints);
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function u = kdeWrapper_x(a,kernelType,probType,bandwidthType,SupportType,nPoints)
      [~,u] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType,NumPoints=nPoints);
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end