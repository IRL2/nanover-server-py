#!/bin/bash
cd protocol
find ./ -name *.proto -print0 | while read -d $'\0' file
do
  echo $file
  python3 -m grpc_tools.protoc -I ./ --python_out=../project/python/narupa-protocol --grpc_python_out=../project/python/narupa-protocol $file
done
