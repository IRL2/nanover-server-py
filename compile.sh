#!/bin/bash
echo "Checking that the requirements.txt is satisfied using pip"
python -m pip install -r ./python-libraries/narupa-core/requirements.txt --user



cd protocol
find ./ -name *.proto -print0 | while read -d $'\0' file


do
  echo $file
  python -m grpc_tools.protoc -I ./ --python_out=../python-libraries/narupa-core --grpc_python_out=../python-libraries/narupa-core $file
done
cd ..
# add init files to all generated python files.
find ./python-libraries/narupa-core/narupa/protocol -type d -exec touch {}/__init__.py \;

cd project/csharp/Narupa.Protocol

dotnet build
