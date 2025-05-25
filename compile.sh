#!/bin/bash

sed -i 's/GPU=1/GPU=0/g' /home/mchao/code/CaLES/build.conf
sed -i 's/NAME := cales/NAME := cales-cpu/g' /home/mchao/code/CaLES/Makefile
make allclean && make libs && make -j


sed -i 's/GPU=0/GPU=1/g' /home/mchao/code/CaLES/build.conf
sed -i 's/NAME := cales-cpu/NAME := cales/g' /home/mchao/code/CaLES/Makefile
make allclean && make libs && make -j
