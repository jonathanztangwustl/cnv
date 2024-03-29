This contains a handful of tools useful for CNV calling (with DNAcopy CBS). Starting with important scripts:
- `indexcov_to_wig.R`: This script converts indexcov output files (indexcov.tar.gz) to WIG files that ichorCNA can input.
- `wig_sizes/offsets`: These are used in the previous script to resize and fix indexcov output to fit ichorCNA input.
- `pon.py`: Used to generate an indexcov panel of normals.
- `run_segs.R`: Used to run DNAcopy segmentation on samples.
- `run_sim_dels.R`: Used to simulate deletions in normal samples (for testing sensitivity/specificity).
- `seg_tools.R`: Tools used for `run_segs.R`.
- `sens_tools.R`: Tools used for `run_sim_dels.R`.
- `sim_del_tools.R`: Tools used for `run_sim_dels.R`.
- `cat_results.R`: Deprecated. For concatenating results from multiple DNAcopy CBS runs.
- `check_id_to_wig.R`: Deprecated. Used to check the conversion from indexcov to WIG, also used to create the WIG offsets.
- 