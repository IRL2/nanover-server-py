
protoc_path="$(pwd)/protoc-3.6.1-win32/bin/protoc.exe"

protocol_directory="$(pwd)/protocol"
csharp_output_directory="$(pwd)/csharp"
cd "protocol"

rm -r "$csharp_output_directory"
mkdir "$csharp_output_directory"

find . -name "*.proto" | while read file
do
  filename="$(pwd)/${file#./}"
  echo "$filename"
  $protoc_path -I="$protocol_directory" --csharp_out="$csharp_output_directory" "$filename"
done


exit 1
