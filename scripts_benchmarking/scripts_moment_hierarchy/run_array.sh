for entry in "./data"/**
do
  ./runopf.jl -f "$entry" >> result.txt
done
