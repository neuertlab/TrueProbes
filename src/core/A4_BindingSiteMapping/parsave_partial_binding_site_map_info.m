function parsave_partial_binding_site_map_info(fname, partial_binding_site_map_info_tmp)
% This function is a wrapper for saving partial_delta_gibson in a parfor loop
  save(fname, 'partial_binding_site_map_info_tmp', '-v7.3')
end