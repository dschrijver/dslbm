#!/bin/bash

TOL=1e-12

h5diff -q -d $TOL ../../data_1.h5 reference/data_1.REF main/u
STATUS_A=$?

h5diff -q -d $TOL ../../data_1.h5 reference/data_1.REF main/v
STATUS_B=$?

h5diff -q -d $TOL ../../data_1.h5 reference/data_1.REF main/phi
STATUS_C=$?

if [ $STATUS_A -ne 0 ] || [ $STATUS_B -ne 0 ] || [ $STATUS_C -ne 0 ]; then
    printf "\e[31mNOT PASSED\e[0m\n"
    touch NOT_PASSED
    exit 1
else
    printf "\e[32mPASS\e[0m\n"
    touch PASS
    exit 0
fi