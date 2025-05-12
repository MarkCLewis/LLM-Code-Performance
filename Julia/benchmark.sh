authors=('ChatGPT4' 'Claude' 'Gemini2.5')

rm times.txt
for author in "${authors[@]}"
do
  cd $author
  for file in *.jl
	do
		echo $author $file >> ../times.txt
		for cnt in {1..7}
		do
			{ time ./run.sh $file ; } 2>> ../times.txt
		done
	done
  cd ..
done