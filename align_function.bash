#!/bin/bash

#basic folder structure for this function: root/raw_data/fw for R1 reads ; root/raw_data/rev for R2 reads;
  #root/bins requires window bedGraphs generated from bedtools makewindows function, named as hg19_chr9_{x}bp.bed, where x= 50 100 200 500 1000 2000 5000; root/submissions for the qsub code files

#expected file: following csv:
   #col1:file name(no 1/2.fq)
   #col2: corrected file name
   #col3: padding+primer length
   #col4: padding sequence
   #col5: start pos of neighbor
   #col6: end pos of neighbor
   #col7: a sequence fragment of undigested product
   #col8-10: R,G,B color for track display
   #col11: use T if allele-specific
   #col12: Name of the allele
   #use comma to separate the names!
#output:
  #root/rev_trimmed:filtered and trimmed reverse end
  #root/fw_matched:forward matched reads
  #root/aligned:aligned bam
  #root/fowbed:bed files converted from bam
    #fowbed/bed:unmodified bed output
    #fowbed/dedup:deduplicated
    #fowbed/no_center:deduped,intrachromosomal,removing reads near the viewpoint 
  #root/bedGraph for binned files
    #bedGraph/bg for background (averaged local)
    #bedGraph/binned for binned with various bin sizes (for reproducibility calculation)
    #bedGraph/window for binned with smaller sliding window. Currently 500bp window/100bp step and 1000bp window/200bp step is used
    #bedGraph/pmf for normalized interaction frequency, generated from 500bp/100bp window file
    #bedGraph/titled for track presentation

#final function: execute indfile rootfolder

#define ref genome desination and target chromosome
refg="/home/ubuntu/ref_genome_masked_snp146"
chrom="chr9"

#takes reverse strand read, filter based on padding sequence, and trim the padding sequence based on chosen length
#will generate tmp intermediate
trimmer(){
    file=$1
    trimstart=$2
    adapt=$3
    awk -v tm=${trimstart} '{print substr($1,1,tm+1)}' $file |awk -v adapt=$adapt '($0~adapt){print NR-1"\n"NR"\n"NR+1"\n"NR+2}' > tmp; awk 'NR==FNR{a[$0]++;next}(FNR in a){print $0}' tmp $file | fastx_trimmer -f $trimstart  
}

matchfwd(){
    revtrimmed=$1
    fw=$2
awk 'NR==FNR&&/^@/{a[$1]++;next;}($1 in a){print FNR"\n"FNR+1"\n"FNR+2"\n"FNR+3}' $revtrimmed $fw > tmp; awk 'NR==FNR {a[$0]++;next;}(FNR in a){print $0}' tmp $fw
}

mergerep(){
    rm *comb*
    for tilt in `ls *_1*`; do
	shortname=`echo $tilt|awk -F "_1" '{print $1}'`
	repnum=`ls ${shortname}*|wc -l`
	for reps in $( seq 1 $repnum ); do
	    alt=`echo $tilt|sed "s/_1/_${reps}/g"`
	    cat $alt >> `echo $tilt|sed 's/_1/_comb/g'`
	done
    done
}

transpose(){
    file=$1
    awk 'BEGIN{OFS="\t";p=0}
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' $file
    }

#function takes a serires of bedGraphs with equal length and combine to a single matrix with first row=title and first column=index
mkbins(){
    chrompath=$1
    outpath=$2
    chromosome=$3
    grep $chromosome $chrompath > tmp.txt
    for i in 50 100 200 500 1000 2000 5000; do
	bedtools makewindows -g tmp.txt -w $i > $outpath/hg19_"$chrom"_${i}bp.bed
	done
	rm tmp.txt
}

bdgtomat(){
    currdir=`pwd`
    folder=$1
    suffix=$2
    output=$3
    cd $folder
    rm tmp*
    rm *comb*txt
    ls *${suffix}|cut -d '.' -f 1 > tmpx
    transpose tmpx > tmp_tit
    for ugly in `ls *${suffix}`  ; do
	awk '{print $4}' $ugly > tmp$ugly
    done
    paste tmp*${suffix} > tmpcomb
    cat tmp_tit tmpcomb > tmpready
    awk 'BEGIN{print 0}{print $2}' $ugly > tmpindex
    paste tmpindex tmpready > $output
    rm tmp*
    cd $currdir
    }

execute(){
    indfile=$1
    rootfolder=$2
    md=`pwd`
    rm -r ${rootfolder}/revtrimmed
    rm -r ${rootfolder}/fw_matched
    rm -r ${rootfolder}/aligned
    rm -r ${rootfolder}/fowbed
    rm -r ${rootfolder}/bedGraph
    mkdir ${rootfolder}/bins
    mkdir ${rootfolder}/revtrimmed
    mkdir ${rootfolder}/fw_matched
    mkdir ${rootfolder}/aligned
    mkdir ${rootfolder}/fowbed
    mkdir ${rootfolder}/fowbed/bed
    mkdir ${rootfolder}/fowbed/dedup
    mkdir ${rootfolder}/fowbed/no_center
    mkdir ${rootfolder}/bedGraph
    mkdir ${rootfolder}/bedGraph/binned
    mkdir ${rootfolder}/bedGraph/window
    mkdir ${rootfolder}/bedGraph/titled
    mkdir ${rootfolder}/bedGraph/bg
    mkdir ${rootfolder}/bedGraph/pmf
    mkdir ${rootfolder}/bedGraph/pmf_norm
    
    rm *tmp*
    rm nr*
    echo name > rc1.tmp
    echo total > rc2.tmp
    echo undigested > rc3.tmp
    echo unique > rc4.tmp

    find ${rootfolder} -name *gz -exec gunzip {} \;

    #create bins
    mkbins ${refg}/hg19_genome_size.txt ${rootfolder}/bins $chrom
    
        for  i in `cat $indfile` ; do
	ori=`echo $i|cut -d ',' -f 1`
	mod=`echo $i|cut -d ',' -f 2`
	len=`echo $i|cut -d ',' -f 3`
	adapt=`echo $i|cut -d ',' -f 4`
	start=`echo $i|cut -d ',' -f 5`
	end=`echo $i|cut -d ',' -f 6`
	udseq=`echo $i|cut -d ',' -f 7`
	color=`echo $i|cut -d ',' -f 8-10`
	AS=`echo $i|cut -d ',' -f 11`
	allele=`echo $i|cut -d ',' -f 12`
	if [ "$AS" = "Y" ]
	then
	    mod1="$mod"_AS_"$allele"
	else
	       mod1="$mod"
	fi

	
#align	
	mv ${rootfolder}/raw_data/rev/${ori}*q ${rootfolder}/raw_data/rev/${mod}_2.fq
	mv ${rootfolder}/raw_data/fw/${ori}*q ${rootfolder}/raw_data/fw/${mod}_1.fq
	
	trimmer ${rootfolder}/raw_data/rev/${mod}_2.fq $len $adapt > ${rootfolder}/revtrimmed/${mod1}_rev.fq
	matchfwd ${rootfolder}/revtrimmed/${mod1}_rev.fq ${rootfolder}/raw_data/fw/${mod}*q > ${rootfolder}/fw_matched/${mod1}_fw.fq
	bwa mem -t 8 ${refg}/hg19_146.fa ${rootfolder}/fw_matched/${mod1}_fw.fq ${rootfolder}/revtrimmed/${mod1}_rev.fq |samtools view -F 4 -Sbh - | samtools sort - ${rootfolder}/aligned/${mod1}
#post align	
	bedtools bamtobed -i ${rootfolder}/aligned/${mod1}.bam > ${rootfolder}/fowbed/bed/${mod1}.bed
	samtools view -h ${rootfolder}/aligned/${mod1}.bam | awk '!a[$3$4]++'|samtools view -F 4 -Sbh -| bedtools bamtobed -i - > tmp_${mod1}.bed
	awk '$5>0&&$6=="+"' tmp_${mod1}.bed | awk '!a[$3]++;!a[$3-1]++{};$a[$3+1]++{}' > tmp1_${mod1}.bed
	awk '$5>0&&$6=="-"' tmp_${mod1}.bed | awk '!a[$4]++;!a[$4-1]++{};$a[$4+1]++{}' > tmp2_${mod1}.bed
	cat tmp1_${mod1}.bed tmp2_${mod1}.bed > ${rootfolder}/fowbed/dedup/${mod1}.bed
	rm tmp*${mod1}.bed
	samtools view -h ${rootfolder}/aligned/${mod1}.bam | awk -v chrom=$chrom -v start=$start -v end=$end '/^@/||($3==chrom&&$4>=end)||($3==chrom&&$4<=start){print $0}' |awk '!a[$3$4]++'|samtools view -F 4 -Sbh -| bedtools bamtobed -i - > tmp_${mod1}.bed
	awk '$5>0&&$6=="+"' tmp_${mod1}.bed | awk '!a[$3]++;!a[$3-1]++{};$a[$3+1]++{}' > tmp1_${mod1}.bed
	awk '$5>0&&$6=="-"' tmp_${mod1}.bed | awk '!a[$4]++;!a[$4-1]++{};$a[$4+1]++{}' > tmp2_${mod1}.bed
	cat tmp1_${mod1}.bed tmp2_${mod1}.bed >${rootfolder}/fowbed/no_center/${mod1}.bed
	rm tmp*${mod1}.bed
	echo $mod1 >> rc1.tmp
	cat ${rootfolder}/revtrimmed/${mod1}* | wc -l >>rc2.tmp
	cat ${rootfolder}/raw_data/rev/${mod}* | grep $udseq | grep $adapt |wc -l >> rc3.tmp
	cat ${rootfolder}/fowbed/no_center/${mod1}*|wc -l >>rc4.tmp
    done
    
paste rc1.tmp rc2.tmp rc3.tmp rc4.tmp > ${rootfolder}/readcount_report.txt

#create comb
    cd ${rootfolder}/fowbed/
    for i in bed dedup no_center; do
	cd ${rootfolder}/fowbed/${i}
	mergerep
	for j in *comb*; do
	    mv $j tmp$j
	    bedtools sort -i tmp$j > $j
	    rm tmp$j
	done
    done
    
#bins
    cd ${rootfolder}/fowbed/no_center
    for i in *bed; do
	for x in 50 100 200 500 1000 2000 5000; do
	    bedf=${rootfolder}/bins/hg19_"$chrom"_${x}bp.bed
	    bedtools intersect -c -e -f 0.50 -F 0.50 -a $bedf -b $i > ${rootfolder}/bedGraph/binned/`echo $i|cut -d '.' -f 1`_${x}bp.bedGraph
	done
    done
#windows
    for i in *bed; do
	for x in 200 450; do
	    bedf=${rootfolder}/bins/hg19_"$chrom"_100bp.bed
	    bedtools window -c -w $x -a $bedf -b $i > ${rootfolder}/bedGraph/window/`echo $i|cut -d '.' -f 1`_`expr ${x} \* 2 + 100`window.bedGraph
	done
    done
$titled
    
#background
    for i in 1200 2450 4950 12450; do
	for j in *bed; do
	    bedtools window -c -w $i -a ${rootfolder}/bins/hg19_"$chrom"_100bp.bed -b $j |awk '{print $4}' > ${rootfolder}/bedGraph/bg/${j}-${i}
	done
    done
    cd  ${rootfolder}/bedGraph/bg
    rm *min
    #based on 500bp min
    for i in `ls|cut -d '-' -f 1|sort -u`; do
	paste -d "\t" ${rootfolder}/bins/hg19_"$chrom"_100bp.bed ${i}-1200 ${i}-2450 ${i}-4950 ${i}-12450|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/5,$5/10,$6/20,$7/50}'|awk 'BEGIN{OFS="\t"}{min=$4;for(i=4;i<=NF;i++)if(min>$i)min=$i;print $1,$2,$3,min}' > ${i}_min;
    done
    rm *.bed-*
    
#pmf*10,000
    cd ${rootfolder}/bedGraph/window
    for i in *500window*ph; do
	awk 'BEGIN{a=0;OFS="\t"}NR==FNR{a+=$4;next}{printf "%s\t%s\t%s\t%f\n",$1,$2,$3,$4/a*10000}' $i $i > ../pmf/pmf_${i};
    done
    cd ../pmf
    rm *txt
    bdgtomat . bedGraph ./pmf_combined.txt
        
#title
    cd ${rootfolder}/bedGraph/pmf
    for  i in *comb*ph ; do
	tit=`echo $i|cut -d '.' -f 1|sed 's/pmf_//g'|sed 's/_comb_500window//g'`
	color=`cat ${md}/${indfile}|grep $tit|cut -d ',' -f 8-10|perl -p -e 's/\r//g'`
	color=`echo $color|awk '{print $1}'`
	lim=`awk 'BEGIN{a=0}(a<$4){a=$4}END{printf "%f",a/3}' $i`
	awk -v tit="$tit" -v c="$color" -v lim="$lim" 'BEGIN{printf "browser position chr9:19480000-21200000\ntrack type=bedGraph name=%s visibility=2 color=%s autoScale=off alwaysZero=on viewLimits=0:%s\n",tit,c,lim}($4!=0){print $0}' $i > ../titled/`echo $i|sed 's/pmf/titled/g'`; done
    cd ../titled
    cat titled* > title_combined.bedGraph
    cd $md
    rm *tmp*
    rm *nr*
    cd ${rootfolder}
#tracks for washu
    cp -r ${rootfolder}/bedGraph/pmf ${rootfolder}/bedGraph/washu
    cd ${rootfolder}/bedGraph/washu
    for i in *.bedGraph; do
	mv $i > tmp$i
	awk '$4>0' tmp$i > $i
	rm tmp $i
	bgzip $i
	tabix -p bed ${i}.gz
    done
    cd ${rootfolder}
    find . -name *fq -exec gzip {} \;
}
    
