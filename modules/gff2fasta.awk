#!/usr/bin/awk -f

BEGIN {
    FS = "\t";
    OFS = "";
    # Define the codon-to-amino acid mappings
    codons["ATA"]="I"; codons["ATC"]="I"; codons["ATT"]="I"; codons["ATG"]="M";
    codons["ACA"]="T"; codons["ACC"]="T"; codons["ACG"]="T"; codons["ACT"]="T";
    codons["AAC"]="N"; codons["AAT"]="N"; codons["AAA"]="K"; codons["AAG"]="K";
    codons["AGC"]="S"; codons["AGT"]="S"; codons["AGA"]="R"; codons["AGG"]="R";
    codons["CTA"]="L"; codons["CTC"]="L"; codons["CTG"]="L"; codons["CTT"]="L";
    codons["CCA"]="P"; codons["CCC"]="P"; codons["CCG"]="P"; codons["CCT"]="P";
    codons["CAC"]="H"; codons["CAT"]="H"; codons["CAA"]="Q"; codons["CAG"]="Q";
    codons["CGA"]="R"; codons["CGC"]="R"; codons["CGG"]="R"; codons["CGT"]="R";
    codons["GTA"]="V"; codons["GTC"]="V"; codons["GTG"]="V"; codons["GTT"]="V";
    codons["GCA"]="A"; codons["GCC"]="A"; codons["GCG"]="A"; codons["GCT"]="A";
    codons["GAC"]="D"; codons["GAT"]="D"; codons["GAA"]="E"; codons["GAG"]="E";
    codons["GGA"]="G"; codons["GGC"]="G"; codons["GGG"]="G"; codons["GGT"]="G";
    codons["TCA"]="S"; codons["TCC"]="S"; codons["TCG"]="S"; codons["TCT"]="S";
    codons["TTC"]="F"; codons["TTT"]="F"; codons["TTA"]="L"; codons["TTG"]="L";
    codons["TAC"]="Y"; codons["TAT"]="Y"; codons["TAA"]="*"; codons["TAG"]="*";
    codons["TGC"]="C"; codons["TGT"]="C"; codons["TGA"]="*"; codons["TGG"]="W";
}

function reverse_complement(seq, rev_comp) {
    n = length(seq);
    rev_comp = "";
    for (i = n; i > 0; i--) {
        base = substr(seq, i, 1);
        if (base == "A") rev_comp = rev_comp "T";
        else if (base == "T") rev_comp = rev_comp "A";
        else if (base == "C") rev_comp = rev_comp "G";
        else if (base == "G") rev_comp = rev_comp "C";
    }
    return rev_comp;
}

function translate(seq, protein) {
    protein = "";
    for (i = 1; i <= length(seq) - 2; i += 3) {
        codon = substr(seq, i, 3);
        protein = protein (codons[codon] ? codons[codon] : "X");
    }
    return protein;
}

/^>/ {
    if (seq_name) sequences[seq_name] = seq;
    seq_name = substr($0, 2);
    seq = "";
    next;
}

!/^>/ {
    seq = seq $0;
}

END {
    if (seq_name) sequences[seq_name] = seq;

    while ((getline < ARGV[2]) > 0) {
        if ($3 == "CDS") {
            split($9, attr, ";");
            id = attr[1]; # Use the attribute as the unique ID for sequences
            if (!cds_seqs[id]) cds_seqs[id] = "";
            cds_seqs[id] = cds_seqs[id] substr(sequences[$1], $4, $5 - $4 + 1);
            strands[id] = $7;
        }
    }

    for (id in cds_seqs) {
        if (strands[id] == "-") {
            cds_seqs[id] = reverse_complement(cds_seqs[id]);
        }
        print ">" id "\n" translate(cds_seqs[id]);
    }
}
