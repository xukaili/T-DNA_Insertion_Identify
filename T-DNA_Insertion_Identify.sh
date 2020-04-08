# put all Fastq files and T-DNA.fa genome.fa in a same dir
# echo 'sh T-DNA_Insertion_Identify.sh  pC1305_T-DNA.fa  genome.fa' |qsub -V -cwd -j y -q bio05.q -o pC.out.txt
# echo 'sh T-DNA_Insertion_Identify.sh  p8_T-DNA.fa      genome.fa' |qsub -V -cwd -j y -q bio05.q -o p8.out.txt

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
done

cat ref.fa  | perl -e '$/ = ">";while(<>){ chomp; my($head,$seq) = split(/\n/,$_,2); next unless($head && $seq); $seq  =~ s/\s+//g; $seq  =~ s/\n*\r*//g; $len=length($seq); print "$head\t1\t$len\n";}' > genome.bed

for i in *.sort.all.bam
do
    samtools index -@ 50  $i
    mkdir  ${i%%.sort.all.bam}
    bamdst -p genome.bed  -o ./${i%%.sort.all.bam} $i
done

sleep 3
grep 'Target] Average depth\t' ./*/coverage.report > depth.txt
rm *cov.wig
# IGV see the T-DNA Insertion Loci
rm   *sort.all.bam*   *cov.wig
