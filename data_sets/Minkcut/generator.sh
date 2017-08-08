#bin/bash

for Nodes in {10..200..10};
do ./rudy -rnd_graph ${Nodes} 10 0 > random10_${Nodes}.txt; done
