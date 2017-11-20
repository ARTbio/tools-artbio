#!/usr/bin/bash


getParams() {
  # Get parameters
  tool_path=$2
  queryString=$4
  dbname=$6
  logfile=$8
  outfile=$11
  params=(tool_path queryString dbname logfile outfile)
  return params
}

params=getParams;
echo $params;
