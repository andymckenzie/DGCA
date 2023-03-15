## Resubmission
This is a resubmission. In this version I have:

* Updated getOption("stringsAsFactors") to getOption("stringsAsFactors", FALSE) to maintain compatability. 

## Explanation of NOTE
* On Windows, I got a NOTE that "Package suggested but not available for checking: 'doMC'". doMC is optionally used for parallelization by one DGCA suggested package (MEGENA), which can be optionally used in DGCA in one function (ddMEGENA), and is not available on Windows.
