source_file=$1
exec_file="${source_file%.*}"
GOMAXPROCS=48 ./$exec_file