#!/bin/bash

#GRENEDALF=/home/lucas/Dropbox/GitHub/grenedalf/bin/grenedalf
GENESIS_APP=/home/lucas/Dropbox/GitHub/genesis/bin/apps/theta_pi_within_pool_sync_jeff

mkdir -p pi
mkdir -p pi_log
rm -f pi/*
rm -f pi_log/*

date

# prepare output table
OUTFILE=sim_results_pi.csv
echo "fname,true_pairwise_het,true_num_seg,pool_size,read_depth,seq_error,num_sites,est_pi_within" > ${OUTFILE}

for line in `tail -n+2 sim_metadata_pi.csv` ; do  
    # split the line
    # fname,true_pairwise_het,true_num_seg,pool_size,read_depth,seq_error,num_sites
    fname=`echo ${line} | cut -f1 -d,`
    true_pairwise_het=`echo ${line} | cut -f2 -d,`
    true_num_seg=`echo ${line} | cut -f3 -d,`
    pool_size=`echo ${line} | cut -f4 -d,`
    read_depth=`echo ${line} | cut -f5 -d,`
    seq_error=`echo ${line} | cut -f6 -d,`
    num_sites=`echo ${line} | cut -f7 -d,`

    # run the computation
    echo "At ${fname%.*}"
    #${GRENEDALF} diversity --sync-path all_sims/${fname} --window-width ${num_sites} --pool-sizes ${pool_size} --min-coverage 3 --measure all --out-dir "pi" --file-prefix "${fname%.*}-" > "pi_log/${fname%.*}.log"
    ${GENESIS_APP} "all_sims/${fname}" ${num_sites} ${pool_size} "${fname%.*}-" > "pi_log/${fname%.*}.log"

    # Get results from grenedalf output file, and attach to big result table
    pi=`cat "pi/${fname%.*}-diversity.csv" | tail -n 1 | cut -f7`
    
    # Get results from grenedalf output file, and attach to big result table
    echo "${line},${pi}" >> ${OUTFILE}

done

date
