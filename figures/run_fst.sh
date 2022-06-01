#!/bin/bash

GRENEDALF=/home/lucas/Dropbox/GitHub/grenedalf/bin/grenedalf

mkdir -p fst
mkdir -p fst_log
rm -f fst/*
rm -f fst_log/*

# prepare output table
OUTFILE=sim_results_fst.csv
echo "fname,true_pi_w,true_pi_b,true_pi_t,true_hudson_fst,true_nei_fst,pool_size,read_depth,seq_error,num_sites,est_spence_nei,est_spence_hudson,est_kofler,est_karlsson" > ${OUTFILE}

for line in `tail -n+2 sim_metadata_fst.csv` ; do  
    # split the line
    # fname,true_pi_w,true_pi_b,true_pi_t,true_hudson_fst,true_nei_fst,pool_size,read_depth,seq_error,num_sites
    fname=`echo ${line} | cut -f1 -d,`
    true_pi_w=`echo ${line} | cut -f2 -d,`
    true_pi_b=`echo ${line} | cut -f3 -d,`
    true_pi_t=`echo ${line} | cut -f4 -d,`
    true_hudson_fst=`echo ${line} | cut -f5 -d,`
    true_nei_fst=`echo ${line} | cut -f6 -d,`
    pool_size=`echo ${line} | cut -f7 -d,`
    read_depth=`echo ${line} | cut -f8 -d,`
    seq_error=`echo ${line} | cut -f9 -d,`
    num_sites=`echo ${line} | cut -f10 -d,`

    # just echo out the line as is, but no new line yet
    echo -n ${line} >> ${OUTFILE}

    # run the computation
    echo "At ${fname%.*} $num_sites"
    for method in spence-nei spence-hudson kofler karlsson ; do
        ${GRENEDALF} fst --sync-path all_sims/${fname} --window-width ${num_sites} --pool-sizes ${pool_size} --method "${method}" --out-dir "fst" --file-prefix "${fname%.*}-${method}-" > "fst_log/${fname%.*}-${method}.log"
        
        # Get results from grenedalf output file, and attach to big result table
        snps=`cat "fst/${fname%.*}-${method}-fst.csv" | tail -n 1 | cut -f4`
        fst=`cat "fst/${fname%.*}-${method}-fst.csv" | tail -n 1 | cut -f5`
        echo -n ",${fst}" >> ${OUTFILE}
    done
    echo >> ${OUTFILE}

done

