mkdir temp_celloracle
for tree in 1 3 5
do
    for ngenes in 110 200 500
    do
        for ncells in 500 800
        do 
            for sigma in 0.1 0.5
            do 
                for seed in 1 2 3 4
                do 
                    data_path=tree${tree}_${ncells}_cells${ngenes}_genes_sigma${sigma}_${seed}
                    rm -rf /project/hli691/scm_paper/1/${data_path}/celloracle_results
                    nohup python -u script_celloracle.py /project/hli691/scm_paper/1/${data_path}/ > temp_celloracle/${data_path}.out


                done  
            done
        done
    done
done