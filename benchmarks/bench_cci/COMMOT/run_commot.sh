for dataset in 1 2 3 4 5 6 7 8
do
    for K in 20 50 100
    do
        for n_spurious in 0 10 30
        do
            echo cci_sc_${dataset}
            echo ${n_spurious} 
            echo ${K}
            python script_commot.py cci_sc_${dataset} ${n_spurious} ${K}
        done
    done 
done