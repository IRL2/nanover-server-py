#!/bin/bash
cd protocol
find ./ -name *.proto -print0 | while read -d $'\0' file
do
  echo $file
  C:/Users/IT036394-Admin/Anaconda3/envs/narupa-protocol/python.exe -m grpc_tools.protoc -I ./ --python_out=../project/python/narupa-protocol --grpc_python_out=../project/python/narupa-protocol $file
done
