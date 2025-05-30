function parsave_partial_on_off(fname, partial_on_off_tmp)
% This function is a wrapper for saving partial_on_off in a parfor loop
  save(fname, 'partial_on_off_tmp', '-v7.3')
end