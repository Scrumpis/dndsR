/* Parametric bootstrap for RELAX (pilot). Loads a RELAX alternative fit saved with --save-fit,
   constrains K = 1 (the RELAX null), refits, reports the null log-likelihood, and writes
   `nsim` alignments simulated from the fitted null model as FASTA.
   placeholders @FIT@, @OUT@, @NSIM@ are filled in per gene */
fit_path = "@FIT@"; out_prefix = "@OUT@"; nsim = @NSIM@;
ExecuteAFile (fit_path);
GetString (lf_name, LikelihoodFunction, 0);
relax.K := 1;
Optimize (res, ^lf_name);
fprintf (stdout, "NULL_LNL\t", Format(res[1][0], 20, 12), "\n");
DATA_FILE_PRINT_FORMAT = 9; /* FASTA */
for (i = 1; i <= nsim; i += 1) {
    DataSet sim = SimulateDataSet (^lf_name);
    DataSetFilter simf = CreateFilter (sim, 1);
    fprintf (out_prefix + i + ".fas", CLEAR_FILE, simf);
}
fprintf (stdout, "DONE\n");
