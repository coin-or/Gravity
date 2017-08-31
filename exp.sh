#bin/bash

for Nodes in {13..15};
#do ./rudy -rnd_graph ${Nodes} 50 0 > random10_${Nodes}.txt; done
#do ../bin/min_k_part ../data_sets/Minkcut/random10_${Nodes}.txt 4 MIP_tree true;
#do ../bin/min_k_part ../data_sets/Minkcut/random10_${Nodes}.txt 3 Node_edge true;
#do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 4 Node_edge false;
#do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 4 MIP_tree false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 SDP false;
#do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 SDP_tree false;
#do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_4.txt 4 MIP_tree false;
#do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_3.txt 3 Node_edge false;
done
