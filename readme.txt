1. Compile
g++ -std=c++11 LORA.cpp -o LORA -O2
g++ -std=c++11 DFS_prune.cpp -o DFS_prune -O2
g++ -std=c++11 HSP.cpp -o HSP -O2

2. Run
./LORA  100000 1 5 5 0 1 k_avg 0.5 1.5 1e-9
./DFS_prune 100000 1 5 5 0 1 0.5 1.5 k_avg
./HSP 100000 1 5 5 0 1 0.5 1.5 k_avg
