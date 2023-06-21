#!/bin/bash

Rscript parrec__genparams.R 2 &
wait
Rscript parrec__genparams.R 4 &