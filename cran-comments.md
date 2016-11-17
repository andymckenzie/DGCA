## Resubmission
This is a resubmission. In this version I have:

* Removed the test that failed on CRAN.

* Enclosed the DOI in the description with angle brackets, i.e. <doi...>.

## Explanation of NOTE
* On win-builder, I got a NOTE that "Package suggested but not available for checking: 'doMC'". doMC is optionally used for parallelization by one DGCA suggested package (MEGENA), which can be optionally used in DGCA in one function (ddMEGENA), and is not available on Windows.
