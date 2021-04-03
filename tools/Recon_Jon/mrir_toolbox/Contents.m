% MRIR_TOOLBOX
%
% Files
%   mrir_array_accelerated_aliasmap                    -
%   mrir_array_accelerated_folding_operator            -
%   mrir_array_birdcage                                -
%   mrir_array_combine                                 -
%   mrir_array_combine_noise_correction_apply          -
%   mrir_array_combine_noise_correction_lookup         -
%   mrir_array_combine_noise_correction_plots          -
%   mrir_array_combine_noise_correction_table          -
%   mrir_array_combine_noise_correction_table_generate -
%   mrir_array_combine_optimalSNR                      -
%   mrir_array_combine_rss                             -
%   mrir_array_GRAPPA                                  - this will be a high-level function that accepts the accelerated data, ACS
%   mrir_array_GRAPPA_1d__plot_dat_acs                 -
%   mrir_array_GRAPPA_1d__plot_kernel                  -
%   mrir_array_GRAPPA_1d_kernel                        -
%   mrir_array_GRAPPA_1d_recon                         -
%   mrir_array_GRAPPA_2d                               -
%   mrir_array_GRAPPA_2d_artifact                      -
%   mrir_array_GRAPPA_2d_fast_griswold                 - wrapper around "fast_grappa3D"
%   mrir_array_GRAPPA_2d_kernel                        -
%   mrir_array_GRAPPA_2d_recon                         -
%   mrir_array_GRAPPA__pick_kernels                    -
%   mrir_array_GRAPPA_acceleration_synthesize          -
%   mrir_array_GRAPPA_check                            -
%   mrir_array_GRAPPA_check_density                    -
%   mrir_array_GRAPPA_gfactor_analytical               -
%   mrir_array_GRAPPA_indices                          -
%   mrir_array_GRAPPA_plottest                         -
%   mrir_array_GRAPPA_prune                            - extract data and ACS lines from sparse arrays
%   mrir_array_SENSE                                   -
%   mrir_array_SENSE_1d_forwardsim                     - simulate sensitivity encoding (for debugging)
%   mrir_array_SENSE_condition_number                  - condition number of SENSE reconstruction
%   mrir_array_SENSE_gfactor_1d                        -
%   mrir_array_SENSE_gfactor_1d__general               -
%   mrir_array_SENSE_gfactor_1d__general_regularize    -
%   mrir_array_SENSE_gfactor_1d__special               -
%   mrir_array_SENSE_gfactor_1d__special_regularize    -
%   mrir_array_SENSE_gfactor_1d_dev                    -
%   mrir_array_SENSE_gfactor_1d_new                    -
%   mrir_array_SENSE_gfactor_2d                        -
%   mrir_array_SENSE_gfactor_2d__general               -
%   mrir_array_SENSE_gfactor_2d__special               -
%   mrir_array_SENSE_gfactor_empirical                 - trim first two time points to remove transients
%   mrir_array_SENSE_gfactor_global_corrcoef           - global correlation coefficient map
%   mrir_array_SENSE_SNR                               -
%   mrir_array_SENSE_Tx_gfactor_1d__general            -
%   mrir_array_SENSE_Tx_gfactor_1d__special            -
%   mrir_array_SNR                                     -
%   mrir_array_SNR_individual                          -
%   mrir_array_SNR_optimal                             -
%   mrir_array_SNR_timeseries                          -
%   mrir_array_stats_matrix                            -
%   mrir_array_stats_sim_data                          -
%   mrir_array_stats_sim_randmtx                       -
%   mrir_array_whitening_apply                         -
%   mrir_array_whitening_operator                      -
%   mrir_artifact_ghost_compute                        -
%   mrir_artifact_ghost_compute__DEBUG                 - "DEBUG" must be set if this function is called. to avoid duplicate
%   mrir_artifact_ghost_compute_TEST                   -
%   mrir_artifact_ghost_correct                        -
%   mrir_artifact_ghost_correct_TEST                   -
%   mrir_artifact_ghost_manual                         -
%   mrir_artifact_phase_stabilization                  -
%   mrir_artifact_t2star_compute                       -
%   mrir_artifact_t2star_plot                          -
%   mrir_conventional                                  - reconstructs conventional (cartesian) acquisitions
%   mrir_conventional_2d                               - reconstructs conventional (cartesian) acquisitions
%   mrir_conventional_3d                               - reconstructs conventional (cartesian) acquisitions
%   mrir_data_zeropad                                  - zero pad arrays to fill in skipped lines/partitions
%   mrir_DIM_AVE                                       -
%   mrir_DIM_CHA                                       -
%   mrir_DIM_COL                                       -
%   mrir_DIM_ECO                                       -
%   mrir_DIM_IDA                                       -
%   mrir_DIM_IDB                                       -
%   mrir_DIM_IDC                                       -
%   mrir_DIM_IDD                                       -
%   mrir_DIM_IDE                                       -
%   mrir_DIM_LIN                                       -
%   mrir_DIM_PAR                                       -
%   mrir_DIM_PHS                                       -
%   mrir_DIM_REP                                       -
%   mrir_DIM_SEG                                       -
%   mrir_DIM_SET                                       -
%   mrir_DIM_SLC                                       -
%   mrir_display_covariance                            - NOT YET WORKING -- colorbar issue
%   mrir_display_gfactor                               -
%   mrir_display_movie                                 -
%   mrir_epi                                           -
%   mrir_epi_birdcage_recon_kawin                      - reconstruct segmented EPI and synthesize 1st birdcage mode
%   mrir_epi_dti                                       -
%   mrir_epi_simple                                    -
%   mrir_evi                                           -
%   mrir_evi_distort_B0_phaseref                       -
%   mrir_evi_GRAPPA                                    -
%   mrir_evi_GRAPPA_ice                                -
%   mrir_evi_GRAPPA_prep                               - perform phase correction and regrid if ramp sampling
%   mrir_fDFT_freqencode                               - forward Discrete Fourier Transform along freq encode
%   mrir_fDFT_phasencode                               - forward Discrete Fourier Transform along phase encode
%   mrir_ice_dimensions                                -
%   mrir_iDFT                                          - inverse Discrete Fourier Transform
%   mrir_iDFT_freqencode                               - inverse Discrete Fourier Transform along freq encode
%   mrir_iDFT_partencode                               - inverse Discrete Fourier Transform along partition encode
%   mrir_iDFT_phasencode                               - inverse Discrete Fourier Transform along phase encode
%   mrir_image_antialias                               - apply anti-aliasing smoothing akin to SYNGO's "IPT"
%   mrir_image_crop                                    - crop image volume reconstructed from oversampled k-space
%   mrir_image_display_aspect                          -
%   mrir_image_display_mag                             -
%   mrir_image_mosaic                                  -
%   mrir_image_slice_deinterleave                      - sort interleaved slices into ascending order
%   mrir_image_vox2ras                                 -
%   mrir_intensity_scale_bin                           -
%   mrir_kspace_trajectory                             -
%   mrir_multishot_segment_collapse                    -
%   mrir_noise_bandwidth                               -
%   mrir_noise_sigma                                   -
%   mrir_partial_echo                                  -
%   mrir_partial_fourier                               - partial fourier reconstruction in phase encoding direction
%   mrir_read_siemens_CxMatrix                         -
%   mrir_read_siemens_simulation                       - loads all "WriteToFile_*.ima" files generated by an ICE simulation
%   mrir_read_siemens_WriteToFile                      -
%   mrir_regrid_trapezoid                              -
%   mrir_regrid_trapezoid_apply                        -
%   mrir_regrid_trapezoid_prep                         -
%   mrir_regrid_trapezoid_rolloff                      -
%   mrir_sensitivity_map                               - estimate coil sensitivity maps from image data
%   mrsr_spectrum_fid                                  -
%   read_meas_prot                                     - read pulse sequence protocol from VB13A-style "meas.dat"
%   read_meas_prot__struct                             -
