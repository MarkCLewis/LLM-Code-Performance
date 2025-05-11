authors=('ChatGPT4' 'Claude' 'Copilot' 'Gemini' 'Gemini2.5')

rm times.txt
for author in "${authors[@]}"
do
  cd $author
  for file in *.c
	do
    ./build.sh $file
		if [ $? -ne 0 ]; then
			echo $author $file >> times.txt
			for cnt in {1..7}
			do
				{ time ./a.out ; } 2>> ../times.txt
			done
		fi
	done
  cd ..
done