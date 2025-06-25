function parsave_probe_target_flanking_info_tmp(fname, probe_target_flanking_info_tmp)
% This function is a wrapper for saving gene_hits_table in a parfor loop
  save(fname, 'probe_target_flanking_info_tmp', '-v7.3')
end