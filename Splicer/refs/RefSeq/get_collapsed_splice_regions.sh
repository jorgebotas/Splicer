
splice_coord_file=$1
outfile=$2
rm $outfile

tmpdir="/tmp/collapsed_splice_regions"
rm -rf $tmpdir
mkdir $tmpdir

genes=$(cut -f 8,8 $splice_coord_file | sort | uniq)
for gene in $genes; do
{
    # Get gene exon/intron boundaries
    grep '\t'$gene'$' $splice_coord_file > $tmpdir/$gene.txt

    wc -l $tmpdir/$gene.txt
    # Merge boundaries and append to outfile
    mergeBed -i $tmpdir/$gene.txt -c 5,7,8 -o distinct >> $outfile

    rm $tmpdir/$gene.txt
}
done

# rm -rf $tmpdir
