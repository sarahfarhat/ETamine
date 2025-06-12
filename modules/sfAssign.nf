#!/usr/bin/env nextflow

process sfAssign {
    tag "$species"
    
    input:
    tuple path(blast_file), path(ltr_fastafile), val(species)


    publishDir '02_LTRAssignmentResults', mode: 'copy'

    output:
    tuple path("${species}.all.copia.fa"), path("${species}.all.gypsy.fa"), path("${species}.all.bel.fa"), val(species)

    script:
    """
    #1: Get the 10 first hits of the blast (keeping only best match per couple of hits). If less than 5 hits remove lines.
    awk '{print \$1}' ${blast_file} | sort -u > ${species}.LTRHDB.blastout6.id
    for i in `cat ${species}.LTRHDB.blastout6.id`; do 
        grep \$i ${blast_file} | sort -k12,12nr | awk '{if(id[\$2]!=1){print \$0; id[\$2]=1}}'; 
    done | awk '{if(id==\$1){if(count<10){print \$0; count=count+1}}else{print \$0; id=\$1; count=1}}' >> ${species}.LTRHDB.blastout6.top10
    awk '{print \$1}' ${species}.LTRHDB.blastout6.top10 | sort | uniq -c | awk '\$1<5 {print \$2}' > ${species}.to_remove
    grep -v -f ${species}.to_remove ${species}.LTRHDB.blastout6.top10 > ${species}.LTRHDB.blastout6.top10.2; mv ${species}.LTRHDB.blastout6.top10.2 ${species}.LTRHDB.blastout6.top10

    #2: Create 2 files. First with ID and all possible assignation and second only the once having more than one assignation
    awk '{split(\$2,t,"_"); print \$1,t[1]}' ${species}.LTRHDB.blastout6.top10 | sort -u > ${species}.filtered.assignation.id
    awk '{split(\$2,t,"_"); print \$1,t[1]}' ${species}.LTRHDB.blastout6.top10 | sort -u | awk '{print \$1}' | sort | uniq -c | awk '\$1!=1 {print \$2}' > ${species}.filtered.assignation.multiSF

    #3: Create 3 files as all.SF.id with ID of the once having exactly one assignation possible 
    awk 'BEGIN{while(getline < "'${species}'.filtered.assignation.multiSF" > 0){id[\$1]=1}} {if(id[\$1]!=1){if(\$2=="Gypsy"){print \$1 >> "'${species}'.all.gypsy.id"}else{if(\$2=="Bel"){print \$1 >> "'${species}'.all.bel.id"}else{if(\$2=="Copia"){print \$1 >> "'${species}'.all.copia.id"}}}}}' ${species}.filtered.assignation.id

    #4: Get blastout6 lines of multiSF copies
    if [ -s ${species}.filtered.assignation.multiSF ]
    then
        grep -f ${species}.filtered.assignation.multiSF ${species}.LTRHDB.blastout6.top10 > ${species}.multiSF.blastout6
        #5: For each multiSF element, count possible SF and remove the once with at least one DIRS_Nimb assignation
        awk '{split(\$2,t,"_"); print \$1,t[1]}' ${species}.multiSF.blastout6 | sort | uniq -c | awk '{if(NR==1){tag[\$3]=\$1; id=\$2}else{if(\$2==id){tag[\$3]=\$1;}else{printf id"\t"; for(i in tag){printf i"\t"tag[i]"\t"} print ""; id=\$2; delete tag; tag[\$3]=\$1}}}END{printf id"\t"; for(i in tag){printf i"\t"tag[i]"\t"} print "";}' | grep -E "Gypsy|Bel|Copia" > ${species}.all_SF.tsv
        grep "DIRS_Nimb" ${species}.multiSF.blastout6 | awk '{print \$1}' | sort -u > ${species}.seqwithDIRNimb; grep -v -f ${species}.seqwithDIRNimb ${species}.all_SF.tsv > ${species}.all_SF.NoNimb.tsv

        #6: Put in ".to_assign.id" assignations with 8 or 9 same SF
        awk '{if(\$3==8 || \$3==9){print \$1,\$2 >> "'${species}'.to_assign.id"}else{if(\$5==9 || \$5==8){print \$1,\$4 >> "'${species}'.to_assign.id"}else{if(\$7==9 || \$7==8){print \$1,\$6 >> "'${species}'.to_assign.id" }}}}' ${species}.all_SF.NoNimb.tsv

        #7: Add new assignation in SF files
        grep "Copia" ${species}.to_assign.id | awk '{print \$1}' >> ${species}.all.copia.id
        grep "Gypsy" ${species}.to_assign.id | awk '{print \$1}' >> ${species}.all.gypsy.id
        grep "Bel" ${species}.to_assign.id | awk '{print \$1}' >> ${species}.all.bel.id
    fi
    #8: Get fasta files
    awk -F '>' 'NR==FNR{ids[\$0]; next} NF>1{f=(\$2 in ids)} f' ${species}.all.copia.id ${ltr_fastafile} | sed -e '/^>/! s/\\(.*\\)/\\U\\1/; s/>_R_/>/' > ${species}.all.copia.fa
    awk -F '>' 'NR==FNR{ids[\$0]; next} NF>1{f=(\$2 in ids)} f' ${species}.all.gypsy.id ${ltr_fastafile} | sed -e '/^>/! s/\\(.*\\)/\\U\\1/; s/>_R_/>/' > ${species}.all.gypsy.fa
    awk -F '>' 'NR==FNR{ids[\$0]; next} NF>1{f=(\$2 in ids)} f' ${species}.all.bel.id ${ltr_fastafile} | sed -e '/^>/! s/\\(.*\\)/\\U\\1/; s/>_R_/>/' > ${species}.all.bel.fa
    """
}
