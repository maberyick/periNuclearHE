function parsave_cellfeat_peri(fname, nucleiCentroids,nucleiCentroids_stro,nucleiCentroids_epi,nucleiCentroids_bund,feat_perinuclei, feat_epi_perinuclei, feat_stro_perinuclei, feat_bund_perinuclei)
  save(fname, 'nucleiCentroids','nucleiCentroids_stro','nucleiCentroids_epi','nucleiCentroids_bund', 'feat_perinuclei', 'feat_epi_perinuclei', 'feat_stro_perinuclei', 'feat_bund_perinuclei');
end