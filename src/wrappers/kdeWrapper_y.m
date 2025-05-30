<<<<<<< HEAD
<<<<<<< HEAD
function u = kdeWrapper_y(a,kernelType,probType,bandwidthType,SupportType)
      [u,~] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType);
=======
function u = kdeWrapper_y(a,kernelType,probType,bandwidthType,SupportType)
      [u,~] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType);
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function u = kdeWrapper_y(a,kernelType,probType,bandwidthType,SupportType)
      [u,~] = kde(double(a),Bandwidth=bandwidthType,Kernel=kernelType,ProbabilityFcn=probType,Support=SupportType);
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end