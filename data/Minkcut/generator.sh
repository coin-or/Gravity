#/usr/local/bin/bash
for Nodes in {50..250..50};
#do ./rudy -rnd_graph ${Nodes} 10 0 > random10_${Nodes}.txt; done
#do ./rudy -spinglass2pm ${Nodes} ${Nodes} 50 0 > spinglass2pm_${Nodes}${Nodes}.txt; done
do 
    echo "size ${Nodes}";
    sudo ./band  ${Nodes} 3;
done
