#!/bin/bash
sequence=$1
thr=$3
scriptdir=$4
spe=$5
SFam=$6

for round in $(seq 1 $2); do 
    echo `date -u` "Running Family for $1 round ${round}"
    dir=${spe}_${SFam}_round_${round}/
    mkdir -p ${dir}/mafft
    mkdir -p ${dir}/cleaned
    #Run usearch and get clusters id
    echo `date -u` "Running usearch for $1 round ${round} " 
    ${scriptdir}/modules/usearch11.0.667_i86linux32 -cluster_fast ${sequence} -id 0.8 -uc ${dir}/results_${spe}_${SFam}.uc -sort length -strand both -clusters ${dir}/${spe}_${SFam}_clust_ > ${dir}/${spe}_${SFam}.outuclust 2>&1
    grep "^C" ${dir}/results_${spe}_${SFam}.uc | awk '$3>1 {print $2}' > ${dir}/${spe}_${SFam}_clustNotsingle.id
    echo `date -u` "usearch for $1 ${round} done." 

    #if nb clusters more than 10 we forget about the singles. There is enough materials with 10 clusters. if less, we run mafft and indel on singles
    nbclust=`wc -l ${dir}/${spe}_${SFam}_clustNotsingle.id | awk '{print $1}'`;
    echo `date -u` "There is ${nbclust} in ${spe} superfamily ${SFam} round ${round}"
    if (( $nbclust<10))
    then
        echo `date -u` "Running ${spe}.all.${SFam}.${round} Take all singles"
        grep "^C" ${dir}/results_${spe}_${SFam}.uc | awk '$3==1 {print $9}' > ${dir}/${spe}_${SFam}_singl.id
        #extract sequences
        awk -F '>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' ${dir}/${spe}_${SFam}_singl.id ${sequence} > ${dir}/${spe}_${SFam}_singl.fa
        #Run mafft on the singles
        echo `date -u` "Running ${spe}_${SFam}.${round} mafft on single"
        s=`ls ${dir}/${spe}_${SFam}_singl.fa`
        mafft --thread $thr --adjustdirection ${s} > ${dir}/mafft/ali.`basename ${s}`
        ali=${dir}/mafft/ali.`basename ${s}`
        python3  ${scriptdir}/modules/indel.py -f ${ali} -p 80 -o ${dir}/cleaned/`basename ${ali} .fa`.indel.fa;
        echo `date -u` "${spe}_${SFam}.${round} mafft and indel on singles done"
    else
        echo `date -u` "Singles in ${spe} superfamily ${SFam} in round ${round} not processed but fasta file available in ${spe}_${SFam}_singl.fa"
        grep "^C" ${dir}/results_${spe}_${SFam}.uc | awk '$3==1 {print $9}' > ${dir}/${spe}_${SFam}_singl.id
        #extract sequences
        awk -F '>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' ${dir}/${spe}_${SFam}_singl.id ${sequence} | awk '/^>'${spe}'_/{sub(/^>'${spe}'_/,">'${spe}'_'${SFam}'_S" ++i "_")} {print}' | awk '{if($1!~/^>/){gsub("-",""); $0=toupper($0)} print $0}' | sed 's/>_R_/>/'>> ${dir}/${spe}_${SFam}_singl.fa
    fi

    #Run mafft and indels on the clusters
    echo `date -u` "Running ${spe}_${SFam}.${round} mafft on clusters"
    for clust in `cat ${dir}/${spe}_${SFam}_clustNotsingle.id`; do
        mafft --thread $thr --adjustdirection ${dir}/${spe}_${SFam}_clust_${clust} > ${dir}/mafft/ali.${clust}.fa
        python3 ${scriptdir}/modules/indel.py -f ${dir}/mafft/ali.${clust}.fa -p 80 -o ${dir}/cleaned/ali.${clust}.indel.fa;
    done
    echo `date -u` "${spe}_${SFam} ${round} mafft and indel on clusters done"
    #Get only one file for all cleaned sequences for next round
    cat ${dir}/cleaned/*.indel.fa | awk '{if($1!~/^>/){gsub("-",""); $0=toupper($0)} print $0}' | sed 's/>_R_/>/' > ${dir}/cleaned/${spe}_${SFam}_cleaned.fa
    if (( $round != $2 ))
    then
        sequence=${dir}/cleaned/${spe}_${SFam}_cleaned.fa
    fi
    echo `date -u` "Family round_${round} done"
done > family.output 2>&1

#Get fasta sequences with a tag on the cluster and the singles sequences
echo `date -u` "Family get fasta sequence file" >> family.output 2>&1
for file in ${dir}/cleaned/ali.*.indel.fa; do
    if [[ "$file" == *singl* ]]; then
        awk '/^>'${spe}'_/{sub(/^>'${spe}'_/,">'${spe}'_'${SFam}'_S" ++i "_")} {print}' $file | awk '{if($1!~/^>/){gsub("-",""); $0=toupper($0)} print $0}' | sed 's/>_R_/>/'
    else
        nb=$(echo "$file" | grep -oP '(?<=ali\.)\d+(?=\.indel\.fa)')
        awk -v nb="$nb" '/^>'${spe}'_/{sub(/^>'${spe}'_/,">'${spe}'_'${SFam}'_C" nb "_")} {print}' $file | awk '{if($1!~/^>/){gsub("-",""); $0=toupper($0)} print $0}' | sed 's/>_R_/>/'
    fi
done
cp ${dir}/${spe}_${SFam}_singl.fa ${spe}_${SFam}_singl.fa