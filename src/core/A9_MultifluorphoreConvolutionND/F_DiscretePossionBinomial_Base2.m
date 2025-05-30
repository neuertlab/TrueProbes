function [Px2D] = F_DiscretePossionBinomial_Base2(Px2D,p)

if (size(p,1) > 0)
    p1 = p(1,1);p2 = p(2,1);    
    Px2D = ImageKernel(Px2D,p1,p2);
    p = p(:,2:end);  
    Px2D = F_DiscretePossionBinomial_Base2(Px2D,p);
end 

end
function Q = ImageKernel(P,p1,p2)
   Ia = [1-p1 p1; 0 0];
   Ib = [[1-p2 p2]' [0 0]'];
  % 1-p1    p1 +  1-p2  0         1-p1/2-p2/2    p1/2       
  %  0      0      p2   0            p2/2         0   normalized
  % 1-p1-p2  p1
  %   p2     0
   Ic = (Ia+Ib)/sum(Ia+Ib,'all');
   Q = conv2(P,Ic); 
end