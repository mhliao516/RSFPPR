eps_arr=(0.5 0.4 0.3 0.2 0.1)
reps_arr=(5e-9 4e-9 3e-9 2e-9 1e-9)

g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o rsf SFMT.c main.cpp -g

./rsf -d dblp_weighted -algo GROUND_TRUTH_BACK -n 10 -alpha 0.2
./rsf -d dblp_weighted -algo GROUND_TRUTH_BACK -n 10 -alpha 0.01

for eps in "${eps_arr[@]}"
do
    ./rsf -d dblp_weighted -algo BACK -n 10 -alpha 0.2 -e ${eps}
done

for reps in "${reps_arr[@]}"
do
    ./rsf -d dblp_weighted -algo RBACK -n 10 -alpha 0.2 -re ${reps}
done

for eps in "${eps_arr[@]}"
do
    ./rsf -d dblp_weighted -algo BACK -n 10 -alpha 0.01 -e ${eps}
done

for reps in "${reps_arr[@]}"
do
    ./rsf -d dblp_weighted -algo RBACK -n 10 -alpha 0.01 -re ${reps}
done

./rsf -d dblp_weighted -algo COMPARE_RESULTS_BACK -n 10 -alpha 0.2
./rsf -d dblp_weighted -algo COMPARE_RESULTS_BACK -n 10 -alpha 0.01