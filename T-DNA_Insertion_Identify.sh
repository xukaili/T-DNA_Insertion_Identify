# sh T-DNA_Insertion_Identify.sh  T-DNA.fa  genome.fa
cat $1  $2 > ref.fa
bwa index  ref.fa

for i in *_R1.fastq.gz
do
    bwa mem -t 50  ref.fa  $i  ${i%%_R1.fastq.gz}_R2.fastq.gz  |samtools view -bh |samtools sort -@ 5 -m 80G -o ${i%%_R1.fastq.gz}.sort.all.bam
done

head=`head -1 $1 | perl -pe 's/\r*\n*//' | perl -pe 's/>//'`
for i in *.sort.all.bam
do
  samtools view  -@ 50  $i | grep $head  >  ${i%%.sort.all.bam}.$head.sam
done

for i in *.sam
do
  igvtools count  -z 5 -w 25 $i   ${i%%.sam}.cov.wig  ref.fa
  cat ${i%%.sam}.cov.wig  | perl -e 'while(<>){print if $_ =~ /chrom/;$_ =~ s/\n*|\r*//g;@line = split(/\t/,$_);next if ($line[1]+$line[2]+$line[3]+$line[4]+$line[5]) < 5 ;print "$_\n";}' > ${i%%.sam}.IGV.txt
  perl -p -i -e '$/ = ""; s/\r*\n(\d+\t)/\t::\t$1/g; s/variab.*span=25\n//g; s/^\d+.*\n//g; s/\t::\t/\n/g;'  ${i%%.sam}.IGV.txt
  igvtools count -w 1 --bases $i  ${i%%.sam}.bases.cov.wig  ref.fa
  cat ${i%%.sam}.bases.cov.wig  | perl -e 'while(<>){print if $_ =~ /chrom/;$_ =~ s/\n*|\r*//g;@line = split(/\t/,$_);$lxk=$line[1]+$line[2]+$line[3]+$line[4]+$line[5];next if ($lxk) < 5 ;print "$line[0]\t$lxk\n";}' > ${i%%.sam}.bases.IGV.txt
  perl -p -i -e '$/ = ""; s/\r*\n(\d+\t)/\t::\t$1/g; s/variab.*span=1\n//g; s/^\d+.*\n//g;'   ${i%%.sam}.bases.IGV.txt
  perl -p -i -e '@a = split(/\t/,$_); for ($i = 2; $i <= $#a; $i = $i + 3){$lxk=$a[$i + 3] - $a[$i]; if ($lxk > 2 and $lxk < 10000) {$a[0]=~/chrom=(\S+)/;print "\nCandiT-DNA\t$1\tStart: $a[$i]\tEnd: $a[$i + 3]\tLength: $lxk\n";}};s/\t/ /g; s/ :: /\t/g;'  ${i%%.sam}.bases.IGV.txt
done
grep 'CandiT-DNA' ./*bases.IGV.txt
grep 'CandiT-DNA' ./*bases.IGV.txt > CandiT-DNA.$head.txt
# Use IGV see the sam and ref.fa files to identify T-DNA Insertion Loci
# grep 'Chr' ./*/*/5.tdna_insertion.bed > TDNAscan.$head.txt


#测序深度分析
#cat ref.fa  | perl -e '$/ = ">";while(<>){ chomp; my($head,$seq) = split(/\n/,$_,2); next unless($head && $seq); $seq  =~ s/\s+//g; $seq  =~ s/\n*\r*//g; $len=length($seq); print "$head\t1\t$len\n";}' > genome.bed
#for i in *.sort.all.bam
#do
#    samtools index -@ 50  $i
#    mkdir  ${i%%.sort.all.bam}
#    bamdst -p genome.bed  -o ./${i%%.sort.all.bam} $i
#done
#sleep 3
#grep 'Target] Average depth\t' ./*/coverage.report > depth.txt
rm   *sort.all.ba*   *wig   igv.log
