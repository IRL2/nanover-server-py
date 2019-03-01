#!/bin/bash
cd protocol
find ./ -name *.proto -print0 | while read -d $'\0' file
do
  echo $file
  python -m grpc_tools.protoc -I ./ --python_out=../project/python/narupa-protocol --grpc_python_out=../project/python/narupa-protocol $file
done
cd ..
# add init files to all generated python files.
find project/python/narupa-protocol/narupa/protocol -type d -exec touch {}/__init__.py \;

