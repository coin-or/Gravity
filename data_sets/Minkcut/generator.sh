#/usr/local/bin/bash
for Nodes in {300..500..100};
#do ./rudy -rnd_graph ${Nodes} 10 0 > random10_${Nodes}.txt; done
#do ./rudy -spinglass2pm ${Nodes} ${Nodes} 50 0 > spinglass2pm_${Nodes}${Nodes}.txt; done
do 
    echo "size ${Nodes}";
    ./band ${Nodes};
done
