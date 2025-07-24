function parsave_partial_designer_stats(fname, partial_designer_stats_tmp)
% This function is a wrapper for saving parsave_partial_designer_stats in a parfor loop
  save(fname, 'partial_designer_stats_tmp', '-v7.3')
end