function parsave_partial_delta_gibson(fname, partial_delta_gibson_tmp)
% This function is a wrapper for saving partial_delta_gibson in a parfor loop
  save(fname, 'partial_delta_gibson_tmp', '-v7.3')
end