#!/bin/bash

TOL=1e-12

h5diff -q -d $TOL ../../data_10.h5 reference/data_10.REF hydro/rho
STATUS_A=$?

h5diff -q -d $TOL ../../data_10.h5 reference/data_10.REF hydro/rho_RED
STATUS_B=$?

h5diff -q -d $TOL ../../data_10.h5 reference/data_10.REF hydro/rho_BLUE
STATUS_C=$?

h5diff -q -d $TOL ../../data_10.h5 reference/data_10.REF hydro/v
STATUS_D=$?

if [ $STATUS_A -ne 0 ] || [ $STATUS_B -ne 0 ] || [ $STATUS_C -ne 0 ] || [ $STATUS_D -ne 0 ]; then
    printf "\e[31mNOT PASSED\e[0m\n"
    touch NOT_PASSED
    exit 1
else
    printf "\e[32mPASS\e[0m\n"
    touch PASS
    exit 0
fi
