
snr_true = 0:0.1:100;


time [snr_bias001, snr_true001, std_bias001, avg_bias001, avg_meas001] = mrir_array_combine_noise_correction_lookup(snr_true, 001);
save mrir_array_combine_noise_correction_table__001 snr_bias001 snr_true001 std_bias001 avg_bias001 avg_meas001

time [snr_bias002, snr_true002, std_bias002, avg_bias002, avg_meas002] = mrir_array_combine_noise_correction_lookup(snr_true, 002);
save mrir_array_combine_noise_correction_table__002 snr_bias002 snr_true002 std_bias002 avg_bias002 avg_meas002

time [snr_bias004, snr_true004, std_bias004, avg_bias004, avg_meas004] = mrir_array_combine_noise_correction_lookup(snr_true, 004);
save mrir_array_combine_noise_correction_table__004 snr_bias004 snr_true004 std_bias004 avg_bias004 avg_meas004

time [snr_bias008, snr_true008, std_bias008, avg_bias008, avg_meas008] = mrir_array_combine_noise_correction_lookup(snr_true, 008);
save mrir_array_combine_noise_correction_table__008 snr_bias008 snr_true008 std_bias008 avg_bias008 avg_meas008

time [snr_bias012, snr_true012, std_bias012, avg_bias012, avg_meas012] = mrir_array_combine_noise_correction_lookup(snr_true, 012);
save mrir_array_combine_noise_correction_table__012 snr_bias012 snr_true012 std_bias012 avg_bias012 avg_meas012

time [snr_bias016, snr_true016, std_bias016, avg_bias016, avg_meas016] = mrir_array_combine_noise_correction_lookup(snr_true, 016);
save mrir_array_combine_noise_correction_table__016 snr_bias016 snr_true016 std_bias016 avg_bias016 avg_meas016

time [snr_bias032, snr_true032, std_bias032, avg_bias032, avg_meas032] = mrir_array_combine_noise_correction_lookup(snr_true, 032);
save mrir_array_combine_noise_correction_table__032 snr_bias032 snr_true032 std_bias032 avg_bias032 avg_meas032

time [snr_bias064, snr_true064, std_bias064, avg_bias064, avg_meas064] = mrir_array_combine_noise_correction_lookup(snr_true, 064);
save mrir_array_combine_noise_correction_table__064 snr_bias064 snr_true064 std_bias064 avg_bias064 avg_meas064

time [snr_bias096, snr_true096, std_bias096, avg_bias096, avg_meas096] = mrir_array_combine_noise_correction_lookup(snr_true, 096);
save mrir_array_combine_noise_correction_table__096 snr_bias096 snr_true096 std_bias096 avg_bias096 avg_meas096

time [snr_bias128, snr_true128, std_bias128, avg_bias128, avg_meas128] = mrir_array_combine_noise_correction_lookup(snr_true, 128);
save mrir_array_combine_noise_correction_table__128 snr_bias128 snr_true128 std_bias128 avg_bias128 avg_meas128

