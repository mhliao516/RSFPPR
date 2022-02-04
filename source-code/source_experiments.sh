eps_arr=(0.5 0.4 0.3 0.2 0.1)

g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o rsf SFMT.c main.cpp -g

mkdir result
mkdir result/dblp_weighted
cd result/dblp_weighted
mkdir gt fora foral foralv speedppr speedl speedlv gtback back backlv rback
cd ../..

./rsf -d dblp_weighted -algo GEN_QUERY -n 100 -alpha 0.2

./rsf -d dblp_weighted -algo GROUND_TRUTH -n 10 -alpha 0.2

for eps in "${eps_arr[@]}"
do
    ./rsf -d dblp_weighted -algo FORA -n 2 -alpha 0.2 -e ${eps}
    ./rsf -d dblp_weighted -algo SPEEDPPR -n 2 -alpha 0.2 -e ${eps}
done

./rsf -d dblp_weighted -algo PLOT_RESULTS -n 2 -alpha 0.2

./rsf -d dblp_weighted -algo GROUND_TRUTH -n 2 -alpha 0.01

for eps in "${eps_arr[@]}"
do
    ./rsf -d dblp_weighted -algo FORA -n 2 -alpha 0.01 -e ${eps}
    ./rsf -d dblp_weighted -algo SPEEDPPR -n 2 -alpha 0.01 -e ${eps}
done

./rsf -d dblp_weighted -algo PLOT_RESULTS -n 2 -alpha 0.01