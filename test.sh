#!/bin/bash

cargo run --release -- resources/test.txt --dna resources/ecoli_m54.txt --truth resources/train.txt > out
diff -y -Z out resources/test_output.txt
