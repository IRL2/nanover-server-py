#!/bin/bash
cd protocol
find ./ -name *.proto -print0 | while read -d $'\0' file
do
  python3 -m grpc_tools.protoc -I ./ --python_out=../build/python --grpc_python_out=../build/python $file
done
