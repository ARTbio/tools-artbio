#!/usr/bin/bash


function getParams() {
  # Get parameters
  tool_path=$2;
  queryString=$4;
  dbname=$6;
  logfile=$8;
  outfile=$11;
  params=(tool_path queryString dbname logfile outfile);
}

function test_getParams() {
  echo "----- getParams test -----";
  getParams;
  if [ ${#params[@]} -ne 5 ]; then
    echo "Didn't pass the test";
    echo "Expected: 5 parameters. Got: ${#params[@]}";
    exit -1;
  fi
}

declare -a params;
test_getParams;
