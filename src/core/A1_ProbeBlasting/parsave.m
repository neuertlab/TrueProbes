function parsave(fname, x,y)
% This function is a wrapper for saving in a parfor loop
  save(fname, 'x', 'y')
end