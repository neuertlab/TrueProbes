function edges = histWrapper_x(a,Method)
      [~,edges] = histcounts(double(a),'binMethod',Method);
end