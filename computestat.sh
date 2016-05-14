#!/bin/bash

max=10
nbmoy=10
rm -f efficacite.dat speedup.dat time.dat
touch efficacite.dat
touch speedup.dat
touch time.dat



#calcul charge total

rm -f Stat/*.dat
mpirun -n 1 ./chppar

chargetot="0"

read -d $'\x04' chargetot < "Stat/stat000.dat" 
echo "chargetot" $chargetot

minchargetot=$chargetot

#calcul sur nbmoy de la chargetot minimun
for j in `seq 2 $nbmoy`

do
    rm -f Stat/*.dat
    mpirun -n 1 ./chppar
    
    read -d $'\x04' chargetot < "Stat/stat000.dat" 
    minchargetot=$(bc <<< "scale=25;0$minchargetot*(0$minchargetot < 0$chargetot)+0$chargetot*(0$minchargetot > 0$chargetot) +0$chargetot*(0$minchargetot == 0$chargetot) ")
    echo "minchargetot" $minchargetot
    
done

echo 1 1 >> efficacite.dat
echo 1 1 >> speedup.dat
echo $i $minchargetot >> time.dat

#calcul de lefficacite,speedup et time pour des nb de processeurs sup a 1

for i in `seq 2 $max`;

do
    efficacite="0.0"
    speedup="0.0"
    totmaxcharge="0.0"

    minmaxcharge=$chargetot    

    #calcul sur nbmoy de la charge max minimun
    for j in `seq 1 $nbmoy`
    
    do
	rm -f Stat/*.dat
	mpirun -n $i ./chppar
	
	maxcharge="0.0"
	for var in `ls Stat`;
	
	do
	    read -d $'\x04' name < "Stat/$var" 
	    echo "name" $name
	    echo "maxcharge" $maxcharge
	    maxcharge=$(bc <<< "scale=25;0$name*(0$name > 0$maxcharge)+ 0$maxcharge*(0$name < 0$maxcharge)+0$name*(0$name == 0$maxcharge) ")
	    echo "newmaxcharge" $maxcharge
	done
	
	minmaxcharge=$(bc <<< "scale=25;0$minmaxcharge*(0$minmaxcharge < 0$maxcharge)+0$maxcharge*(0$minmaxcharge > 0$maxcharge) +0$maxcharge*(0$minmaxcharge == 0$maxcharge)  ")

    done

    efficacite=$(bc <<< "scale=25;($chargetot/($minmaxcharge*$i))")
    speedup=$(bc <<< "scale=25;($chargetot/($minmaxcharge))")
    time=$(bc <<< "scale=25;$minmaxcharge")

    echo $i $efficacite >> efficacite.dat
    echo $i $speedup >> speedup.dat
    echo $i $time >> time.dat
    

    
    
done
rm -f plotstat.gnu
touch plotstat.gnu

echo "set terminal postscript eps enhanced color" >> plotstat.gnu
echo "set output 'plotefficacite.eps' " >> plotstat.gnu
echo "plot 'efficacite.dat' with linespoint ls 7 " >> plotstat.gnu
echo " " >> plotstat.gnu
echo "set terminal postscript eps enhanced color" >> plotstat.gnu
echo "set output 'plotspeedup.eps' " >> plotstat.gnu
echo "plot 'speedup.dat' with linespoint ls 7, x title 'ideal'" >> plotstat.gnu
echo " " >> plotstat.gnu
echo "set terminal postscript eps enhanced color" >> plotstat.gnu
echo "set output 'plottime.eps' " >> plotstat.gnu
echo "plot 'time.dat' with linespoint ls 7" >> plotstat.gnu
echo " " >> plotstat.gnu

gnuplot plotstat.gnu
