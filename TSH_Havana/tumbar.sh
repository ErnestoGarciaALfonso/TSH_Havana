#!/bin/bash


for i in `qstat | awk '{print $1}'`; do qdel $i; done
