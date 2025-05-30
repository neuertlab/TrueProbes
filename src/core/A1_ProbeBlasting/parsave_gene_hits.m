function parsave_gene_hits(fname, gene_hits_tmp)
% This function is a wrapper for saving gene_hits in a parfor loop
  save(fname, 'gene_hits_tmp', '-v7.3')
end