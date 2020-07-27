# adjustcarryforward

An experiment in re-assessing height measurements initially marked as
`Exclude-Carried-Forward` based on growth velocity is included in the current version of
this package.

Because of the complexity involved in this approach, this is implemented as an
independent function with its own driver script for further study. Our hope is to revise
this strategy, and, if feasible, incorporate it into the main `cleangrowth()` algorithm.

## Running the experiment

The primary function is called `adjustcarryforward()`, and is implemented in the file
`R/adjustcarryforward.R`. It is separate from the main `R/growth.R` file as it is set up
to be run separately, using the output from `cleangrowth()`. The `adjustcarryforward()`
function is available in the main package namespace.

The primary tool from running this function is the script `exec/testadjustcf.R`. This
will run the experimental function on an existing dataset (specified as a CSV file) and
pass it a sweep of parameter values for testing the new strategy. This tool will
produce a combined output file that includes the original result alongside the resulting
re-inclusion determination for each measurement for each parameter value run.

For example, if you have run `cleangrowth()` on the `syngrowth` synthetic data set
provided with `growthcleanr` as described in the main `README.md` file, save it to a CSV
file:

```R
> fwrite(cleaned_data, "cleaned.csv", row.names = F)
```

The sweep script is executed from the command line on the cleaned data file:

```bash
% Rscript exec/testadjustcf.R cleaned.csv
```

By default, the script will generate a range of values with nine steps for the
following parameters, where the min and max surround the default value:

| parameter | default | min | max |
| - | - | - | - |
`minfactor` | 0.5 | 0 | 1
`maxfactor` | 2 | 0 | 4
`banddiff` | 3 | 0 | 6
`banddiff_plus` | 5.5 | 0 | 11
`min_ht.exp_under` | 2 | 0 | 4
`min_ht.exp_over` | 0 | -1 | 1
`max_ht_exp_under` | 0.33 | 0 | 0.66
`max_ht.exp_over` | 1.5 | 0 | 3

The default number of sweep steps is 9; this can be changed with the option
`--gridlength`. For example, for a 9-step sweep, the parameters passed to the
function in each pass will be:

```R
run  minfactor  maxfactor  banddiff  banddiff_plus  min_ht.exp_under  min_ht.exp_over  max_ht.exp_under  max_ht.exp_over
1    0          0          0         0              0                 -1               0                 0
2    0.125      0.5        0.75      1.375          0.5               -0.75            0.0825            0.375
3    0.25       1          1.5       2.75           1                 -0.5             0.165             0.75
4    0.375      1.5        2.25      4.125          1.5               -0.25            0.2475            1.125
5    0.5        2          3         5.5            2                 0                0.33              1.5
6    0.625      2.5        3.75      6.875          2.5               0.25             0.4125            1.875
7    0.75       3          4.5       8.25           3                 0.5              0.495             2.25
8    0.875      3.5        5.25      9.625          3.5               0.75             0.5775            2.625
9    1          4          6         11             4                 1                0.66              3
```

The output, in the `output/` directory will contain the sweep parameters, like
the above, in a file called `params.csv`, and the output with adjustment results
in a file called `all-adjusted.csv`.

For example, a 5-step sweep would be run with this command:

```bash
% Rscript exec/textadjustcf.R --gridlength 5 cleaned.csv
```

The parameter set for the sweep in file `output/params.csv` would be:

```R
run  minfactor  maxfactor  banddiff  banddiff_plus  min_ht.exp_under  min_ht.exp_over  max_ht.exp_under  max_ht.exp_over
1    0          0          0         0              0                 -1               0                 0
2    0.25       1          1.5       2.75           1                 -0.5             0.165             0.75
3    0.5        2          3         5.5            2                 0                0.33              1.5
4    0.75       3          4.5       8.25           3                 0.5              0.495             2.25
5    1          4          6         11             4                 1                0.66              3
```

Note that an odd-numbered length will include the default values in the middle run of the
sweep (hence the examples w/5 and 9 step sweeps).

And the first few result rows in `output/all-adjusted.csv` would be:

```R
id     subjid    sex  agedays  param     measurement  exclude                    run-1      run-2      run-3      run-4      run-5
1510   775155    0    889      HEIGHTCM  84.9         Exclude-Duplicate          Missing    Missing    Missing    Missing    Missing
1511   775155    0    889      HEIGHTCM  89.06        Include                    No Change  No Change  No Change  No Change  No Change
1512   775155    0    1071     HEIGHTCM  92.5         Include                    No Change  No Change  No Change  No Change  No Change
1513   775155    0    1253     HEIGHTCM  96.2         Include                    No Change  No Change  No Change  No Change  No Change
1514   775155    0    1435     HEIGHTCM  96.2         Exclude-Carried-Forward    No Change  No Change  Include    Include    Include
1515   775155    0    1435     HEIGHTCM  99.692       Include                    No Change  No Change  No Change  No Change  No Change
1516   775155    0    1806     HEIGHTCM  106.1        Include                    No Change  No Change  No Change  No Change  No Change
1517   775155    0    2177     HEIGHTCM  112.3        Include                    No Change  No Change  No Change  No Change  No Change
1518   775155    0    889      WEIGHTKG  13.1         Include                    No Change  No Change  No Change  No Change  No Change
```

The fifth row in the example above demonstrates the results of the experimental script;
for runs 1 and 2, the result is not changes, but for runs 3-5, the measurement is
adjusted for reinclusion. To demonstrate the range, the following is an extract of
measurements only marked as carried forward exclusions by `cleangrowth()`:

```R
id     subjid     sex  agedays  param     measurement  exclude                  run-1      run-2      run-3      run-4      run-5
1514   775155     0    1435     HEIGHTCM  96.2         Exclude-Carried-Forward  No Change  No Change  Include    Include    Include
1521   775155     0    1435     WEIGHTKG  15.3         Exclude-Carried-Forward  No Change  No Change  No Change  No Change  No Change
7952   1340377    1    1806     HEIGHTCM  107.1        Exclude-Carried-Forward  No Change  Include    Include    Include    Include
7967   1340377    1    1806     WEIGHTKG  18.4         Exclude-Carried-Forward  No Change  No Change  No Change  No Change  No Change
41775  3643526    1    1253     HEIGHTCM  87.808       Exclude-Carried-Forward  Include    Include    Include    Include    Include
44901  3706097    0    4032     HEIGHTCM  138.8        Exclude-Carried-Forward  No Change  Include    Include    Include    Include
30011  5792371    1    3661     HEIGHTCM  145.4        Exclude-Carried-Forward  No Change  Include    Include    Include    Include
30013  5792371    1    4032     HEIGHTCM  145.4        Exclude-Carried-Forward  No Change  No Change  No Change  No Change  No Change
30016  5792371    1    1071     WEIGHTKG  15.9         Exclude-Carried-Forward  No Change  No Change  No Change  No Change  No Change
```

Some of these values are not adjusted at all; one is from run 1 on, a few are from run 2
on, and one is from run 3 on.
