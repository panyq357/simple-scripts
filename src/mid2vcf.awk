# Convert RiceHap3 mid file to VCF
# Usage: awk -f mid2vcf.awk chr1.mid > chr1.vcf

NR == 1 {
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
    for (i=9;i<NF;i++) {
        printf("%s\t", $i)
    }
    printf("%s\n", $NF)
}

NR > 1 {
    chrom = $1
    pos = $2
    ID = sprintf("%s-%s", chrom, pos)
    ref = $3
    alt = $4
    qual = "."
    filter = "PASS"
    info = "."
    format = "GT"
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
           chrom, pos, ID, ref, alt, qual, filter, info, format)
    for (i=9;i<NF;i++) {
        gt = $i
        if (gt == ref) printf("0/0\t")
        else if (gt == alt) printf("1/1\t")
        else printf("./.\t")
    }
    if ($NF == ref) printf("0/0\n")
    else if ($NF == alt) printf("1/1\n")
    else printf("./.\n")
}

