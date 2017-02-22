#!/bin/sh

##State the custom file working on
echo -e "\t\tWorking on $1"

##Run nextflow pipeline
echo -e "\n =>  Processing your input files ... Sit back & relax ...\n"

#1. VARIANT PIPE
vari=$( echo $0 | sed s/nextVAP.sh/'src\/nextflow\/varipipe.nf'/g )

nextflow run $vari -c $1 -resume > .information
cat .information


#Determine VARIANT PIPE status
status=( $(grep "Execution status" .information | wc -m) )

if [-n "$2"]
  then
  if [$2 == "clean"]
    then
    if [ $status = 26 ]
    then
      echo -e "\tRemoving temporary files"
      rm -rf work/ .information .nextflow.* trace_pipeline.txt*
    fi
  fi
fi
echo
echo -e "\t\t...The End..."


