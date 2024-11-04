#!/bin/bash

## interproscan runs a number of the analyses in parallel, so we should not
## use too many threads or we will use up too many resources
## This may use something like 20 threads or so
ls *.fa | xargs -P 10 -I{} ../run_scans.sh {}
