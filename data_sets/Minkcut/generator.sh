#bin/bash

for Nodes in {10..200};
do ./rudy -rnd_graph ${Nodes} 50 0 > random_${Nodes}.txt; done
