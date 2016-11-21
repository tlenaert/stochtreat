#!/bin/bash

./stochtreat  run=$id resistance=$resistance patients=$patients \
    output=$output treattest=$treattest treattime=$treattime \
    > $id.dat
#epsc=$epsc epsb=$epsb epsn=$epsn 
