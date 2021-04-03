

opt.ReturnStruct = 1;
opt.PhascorCollapseSegments = 0;

meas = read_meas_dat('meas.dat', opt);


evi_raw = meas.patrefscan;
phascor1d = meas.patrefscan_phascor;


evi_raw_apod = mrir_filter_raw_apodize_1d(evi_raw, mrir_DIM_COL, 0.15);
evi_hyb_roft = mrir_iDFT_freqencode(evi_raw_apod);

linear_fit_coeff = mrir_artifact_ghost_compute(phascor1d);
evi_hyb_corr = mrir_artifact_ghost_correct(evi_hyb_roft, linear_fit_coeff);

evi_hyb_coll = mrir_multishot_segment_collapse(evi_hyb_corr, meas.evp);

prot_trapezoid = mrir_regrid_trapezoid_prep(meas.prot, size(evi_raw, 1));
evi_hyb_grid = mrir_regrid_trapezoid(evi_hyb_coll, meas.prot, prot_trapezoid, 1);


evi_hyb_apod = mrir_filter_raw_apodize_1d(evi_hyb_grid, mrir_DIM_LIN, 0.15);
evi_hyb_peft = mrir_iDFT_phasencode(evi_hyb_apod);

evi_hyb_paft = mrir_iDFT_partencode(evi_hyb_peft);

evi_img = mrir_image_crop(evi_hyb_paft);
evi_rss = mrir_array_combine_rss(evi_img);

