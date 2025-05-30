function parsave_gene_hits_table(fname, gene_hits_table_tmp)
% This function is a wrapper for saving gene_hits_table in a parfor loop
  save(fname, 'gene_hits_table_tmp', '-v7.3')
end