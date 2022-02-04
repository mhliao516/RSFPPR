g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o rsf SFMT.c main.cpp -g

./rsf -d dblp_weighted -algo GEN_QUERY -n 100 -alpha 0.2