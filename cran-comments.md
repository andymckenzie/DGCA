## Resubmission
This is a resubmission. In this version I have:

* Removed the class() function calls so that the package will be compatible with R 4.0.0 when it is released. 

## Explanation of NOTE
* On win-builder, I got a NOTE that "Package suggested but not available for checking: 'doMC'". doMC is optionally used for parallelization by one DGCA suggested package (MEGENA), which can be optionally used in DGCA in one function (ddMEGENA), and is not available on Windows.
