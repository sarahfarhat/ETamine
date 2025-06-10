#!/usr/bin/env nextflow

process blast2GFF {
    tag "$species"

    input:
    tuple path(blastout), val(species)
    
    output:
    tuple path("${blastout}.gff"), val(species)

    script :
    """
    awk 'BEGIN{deb[1]=1; fin[1]=1}\
    {OFS="\t"; \
            if(\$1!=id){\
                            delete deb; delete fin;\
                            split(\$2,bm,"-");
                    if(\$7<\$8){\
                            print \$1,"blastRT","CDS",\$7,\$8,".","+",".",bm[1]"_"\$1;\
                            deb[\$1".1"]=\$7; \
                            fin[\$1".1"]=\$8; \
                    }else{\
                            print \$1,"blastRT","CDS",\$8,\$7,".","-",".",bm[1]"_"\$1; \
                            deb[\$1".1"]=\$8; \
                            fin[\$1".1"]=\$7} \
                    id=\$1; \
                    cpt=1; \
            }else{\
                    if(\$7<\$8){\
                            topr="yes"; \
                            for(k in deb){\
                                    if((\$7>=deb[k] && \$8<=fin[k])||(\$7<=deb[k] && \$8>=fin[k])) {topr="no"}\
                                    debdiff=\$7-deb[k]; findiff=\$8-fin[k];\
                                    if(debdiff<0){if((\$8-deb[k])>100) topr="no"}else{if((fin[k]-\$7)>100){topr="no"}} \
                            } \
                            if(topr=="yes"){cpt=cpt+1; print \$1,"blastRT","CDS",\$7,\$8,".","+",".",bm[1]"_"\$1; deb[\$1"."cpt]=\$7; fin[\$1"."cpt]=\$8;}\
                    }else{ \
                            topr="yes"; \
                            for(k in deb){ \
                                    if((\$8>=deb[k] && \$7<=fin[k])||(\$8<=deb[k] && \$7>=fin[k])){topr="no"}\
                                    debdiff=\$8-deb[k]; findiff=\$7-fin[k];\
                                    if(debdiff<0){if((\$8-deb[k])>100) topr="no"}else{if((fin[k]-\$7)>100){topr="no"}} \
                            } \
                            if(topr=="yes"){cpt=cpt+1; print \$1,"blastRT","CDS",\$8,\$7,".","-",".",bm[1]"_"\$1; deb[\$1"."cpt]=\$8; fin[\$1"."cpt]=\$7;}\
                    }\
            }\
    }' ${blastout} | sort -k1,1 -k4,4n > ${blastout}.gff
    """
}

