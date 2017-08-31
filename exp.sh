#bin/bash


# Table 1
for Nodes in {3..16};
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 MIP true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 MIP false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 MIP_tree false;
done

for Nodes in {3..16};
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 MIP true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 MIP false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 MIP_tree false;
done







# Table 2
for Nodes in {50..250..50};
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_3.txt 3 Node_edge true;
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_3.txt 3 Node_edge false;
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_3.txt 3 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_3.txt 3 MIP_tree false;
done

for Nodes in {50..250..50};
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_4.txt 4 Node_edge true;
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_4.txt 4 Node_edge false;
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_4.txt 4 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/band${Nodes}_4.txt 4 MIP_tree false;
done

for Nodes in {10..15};
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 Node_edge true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 Node_edge false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 MIP_tree false;
done

for Nodes in {10..15};
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 4 Node_edge true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 4 Node_edge false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 4 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 4 MIP_tree false;
done

# Table 3
for Nodes in {11..15};
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 SDP true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 SDP false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 SDP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2g_${Nodes}${Nodes}.txt 3 SDP_tree false;
done

for Nodes in {11..15};
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 SDP true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 SDP false;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 SDP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/spinglass2pm_${Nodes}${Nodes}.txt 3 SDP_tree false;
end


# Table 4
for Nodes in {10..150..10};
do ../bin/min_k_part ../data_sets/Minkcut/random10_${Nodes}.txt 3 MIP_tree true;
do ../bin/min_k_part ../data_sets/Minkcut/random10_${Nodes}.txt 3 SDP_tree true;
done
