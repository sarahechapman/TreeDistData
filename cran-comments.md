## Test environments
* Windows 10 on local machine, R 4.0.2
* Windows 10 via `check_win_devel()`, R devel
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0 and release
* Using `check_for_cran()`

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:

> New submission. 

This package is a new submission.


> Possibly mis-spelled words in DESCRIPTION:
>   unrooted (10:73)

False positive: the spelling is correct.

> Availability using Additional_repositories specification:
>   ?   ?   https://ms609.github.io/packages

Repository includes 'TreeTools' v1.1.0, 'TreeDist' v1.1.0 & 'phangorn' v2.6.0.
Latest versions of 'TreeTools' and 'TreeDist' not yet available via all CRAN 
mirrors, but necessary for testing on e.g. `check_for_cran()`.
`phangorn` v2.6.0 recommended for users (as 2.5.5 has memory leaks), but not yet
available through CRAN.

> Size of tarball: 18827486 bytes
> [...]
> * checking installed package size ... NOTE
>   installed size is 34.2Mb
>   sub-directories of 1Mb or more:
>     data  34.0Mb
      
This package contains large datasets used to build the vignettes of
the package 'TreeDist'.  I have compressed the datasets as far
as possible.
  
These datasets should very seldom change, and are primarily used to build 
vignettes (though there are some other use cases);
as such, I feel that it is most appropriate to serve them as a standalone
package.


## Downstream dependencies
There are currently no downstream dependencies for this package.
