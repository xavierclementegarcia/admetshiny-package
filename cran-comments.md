## R CMD check results

0 errors | 0 warnings | 0 note

* Note: New submission 

## Resubmission

Dear Konstanze Lauseker,

This is a resubmission of the admetshiny package. The following
issues raised in the previous review have been addressed:

1. Removed single quotes around acronyms and expanded all acronyms
   in the Description field.
2. Replaced \dontrun{} with \donttest{} or if(interactive()) where
   appropriate; unwrapped examples that execute in <5 sec.
3. Removed set.seed(42) from the plotUMAP() function body.
4. Rewrote the generateReport() example to be self-contained.

The example for getSmilesFromIdentifiers() remains in dontrun{}
because it queries the PubChem REST API.

Duration: 2m 9.8s
R CMD check --as-cran passes with 0 errors, 0 warnings, 0 notes

Best regards,
Xavier Clemente Garcia Cevallos
