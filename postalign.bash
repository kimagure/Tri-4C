#!/bin/bash

#transpose file
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

#three inputs; bdg file name, output filename, bp. To combine usage with combmacs2() use fixed_ for output prefix
fixmacs2(){
    file=$1
    out=$2
    bp=$3
    awk -v bp=$bp 'BEGIN{a=0;OFS="\t"}NR==1{a=0;print $1,a,a+bp,$4;c=$4;a=a+bp;next}(a<$2){while(a!=$2){print $1,a,a+bp,c;a=a+bp}}(a==$2){while(a<$3){b=$3;if(b>a+bp){b=a+bp};print $1,a,b,$4;c=$4;a=a+bp}}' $file | awk '$2>=0{print $0}' > $out
}

#input=folder path. Integrates all fixed_*bdg in to matrix
combmacs2(){
    folder=$1
    output=$2
    cd $folder
    for i in fixed*bdg; do
	awk '{print $4}' $i > tmp$i
    done
    paste tmpfixed* > tmpcomb
    ls fixed*bdg|cut -d '.' -f 1 > tmpx
    transpose tmpx > tmp_tit
    cat tmp_tit tmpcomb > tmpready
    awk 'BEGIN{print 0}{print $2}' $i > tmpindex
    paste tmpindex tmpready > ../$output
    rm tmp*
    cd ..
}


#convert a series of bedgraoh to a matrix
bdgtomat(){
    currdir=`pwd`
    folder=$1
    suffix=$2
    output=$3
    cd $folder
    rm tmp*
    rm *txt
    ls *${suffix}|cut -d '.' -f 1 > tmpx
    transpose tmpx > tmp_tit
    for ugly in `ls *${suffix}`  ; do
	awk '{print $0}' $ugly > tmp$ugly
    done
    paste tmp*${suffix} > tmpcomb
    cat tmp_tit tmpcomb > tmpready
    awk 'BEGIN{print 0}{print $2}' $ugly > tmpindex
    paste tmpindex tmpready > $output
    rm tmp*
    cd $currdir
}



analysis(){

    rootfolder=$1
    cd $rootfolder
    rm -r ${rootfolder}/macs2
    mkdir ${rootfolder}/macs2
    mkdir ${rootfolder}/macs2/ppois
    mkdir ${rootfolder}/macs2/logFE
    mkdir ${rootfolder}/macs2/IDR
    mkdir ${rootfolder}/macs2/IDR/init
    mkdir ${rootfolder}/macs2/IDR/res
    cd ${rootfolder}/bedGraph/window
    for i in *500window*ph; do
	j=`echo $i|sed 's/_500window.bedGraph/.bed_min/g'`
	macs2 bdgcmp -m ppois -t $i -c ../bg/$j --o-prefix ${rootfolder}/macs2/ppois/ppois_`echo $j|cut -d '.' -f 1`
	macs2 bdgcmp -p 5 -m logFE -t $i -c ../bg/$j --o-prefix ${rootfolder}/macs2/logFE/logfe_`echo $j|cut -d '.' -f 1`
    done
    cd ${rootfolder}/macs2/ppois
    rm fixed*
    for i in *bdg; do fixmacs2 $i fixed_${i} 100;done
    cd ${rootfolder}/macs2/logFE
    rm fixed*
    for i in *bdg; do fixmacs2 $i fixed_${i} 100;done
    cd ${rootfolder}/macs2
    combmacs2 logFE logFE_comb.txt
    combmacs2 ppois ppois_comb.txt
    #pk call
    cd ${rootfolder}/macs2/ppois
    for i in fixed*_1_*bdg; do
	j=`echo $i| sed 's/_1_/_2_/g'`
	paste $i $j | awk 'BEGIN{OFS="\t"}$4>0&&$8>0{print $1,$2,$3,NR,$4,".",$4,$4,$4,50}' > ../IDR/init/${i}.bed
	paste $i $j | awk 'BEGIN{OFS="\t"}$4>0&&$8>0{print $1,$2,$3,NR,$8,".",$8,$8,$8,50}' > ../IDR/init/${j}.bed
    done
    cd ${rootfolder}/macs2/IDR/init
    for i in *_1_*bed; do j=`echo $i|sed 's/_1_/_2_/g'`; idr --input-file-type narrowPeak --samples $i $j --initial-mu 1.5 --initial-sigma 0.3 --initial-rho 0.8 --output-file ../res/`echo ${i}|cut -d '.' -f 1|sed 's/_1_/_/g'`_output.peak; done
    cd ${rootfolder}/macs2/IDR/res
    for i in *output.peak; do cat $i|awk '$5>=540'|cut -f 1-5|bedtools sort -i - |bedtools merge -i - > ${i}_merged; done
    #washu track
    for i in *merged; do awk '{print $0"\t"1}' $i > ${rootfolder}/bedGraph/washu/$i ; done
    cd ${rootfolder}/bedGraph/washu
    for i in *_merged; do
	bgzip $i
		tabix -p bed ${i}.gz
    done
    cd ${rootfolder}
}

