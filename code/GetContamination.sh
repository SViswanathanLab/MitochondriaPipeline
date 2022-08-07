#!/bin/bash

echo $1;
echo $2;
grep "SampleID" $1 > $2;
