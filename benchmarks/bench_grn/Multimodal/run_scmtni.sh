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
                    indir=tree${tree}_${ncells}_cells${ngenes}_genes_sigma${sigma}_${seed}
                    filelist=${indir}/filelist.txt
                    regfile=${indir}/regulators.txt

                    # preprocess the count matrix and prior GRN
                    python preprocess_scmtni.py ${indir}
                    # running scmtni
                    python Scripts/PreparescMTNIinputfiles.py --filelist $filelist --regfile $regfile --indir $indir --outdir results_${indir} --splitgene 50 --motifs 1
                    Code/scMTNI -f ${indir}/testdata_config.txt -x50 -l ${indir}/TFs_OGs.txt -n ${indir}/AllGenes.txt -d ${indir}/celltype_tree_ancestor.txt -m ${indir}/testdata_ogids.txt -s ${indir}/celltype_order.txt -p 0.2 -c yes -b -0.9 -q 2 
                    # postprocessing the cluster-specific GRNs, using averaging result
                    python postprocess_scmtni.py ${indir}
                    # save intermid files
                    mv results_${indir} /project/hli691/scm_paper/1/${indir}/scmtni_results/
                    rm -rf ${indir}

                done  
            done
        done
    done
done