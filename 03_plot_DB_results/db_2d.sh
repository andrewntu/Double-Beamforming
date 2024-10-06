#!/bin/bash
gmt gmtset PAPER_MEDIA=a3
gmt gmtset MAP_FRAME_PEN=0.5p
sta=$1
period=$2
if [ $# != 2 ];
then
echo Error Input Format!!
echo 'Example: sh db_2d M02_M26 3'
exit
fi

dist_min=`sort -n -k2 Output/period_${period}/${sta}_${period}.txt | head -n 1 | gawk '{print$2-2.5}'`
dist_max=`sort -n -k2 Output/period_${period}/${sta}_${period}.txt | tail -n 1 | gawk '{print$2+2.5}'`
fund_min_coarse=`sort -n -k1 Output/period_${period}/1st_${sta}_coarse | head -n 1 | gawk '{print$1}'`
fund_max_coarse=`sort -n -k1 Output/period_${period}/1st_${sta}_coarse | tail -n 1 | gawk '{print$1}'`
fund_source_slow_min=`sort -n -k1 Output/period_${period}/1st_${sta}_coarse | head -n 1 | gawk '{print$1-0.1}'`
fund_source_slow_max=`sort -n -k1 Output/period_${period}/1st_${sta}_coarse | tail -n 1 | gawk '{print$1+0.1}'`
fund_rec_slow_min=`sort -n -k2 Output/period_${period}/1st_${sta}_coarse | head -n 1 | gawk '{print$2-0.1}'`
fund_rec_slow_max=`sort -n -k2 Output/period_${period}/1st_${sta}_coarse | tail -n 1 | gawk '{print$2+0.1}'`
if [ "`ls Output/period_${period}/2nd_${sta}_coarse | wc -l`" -ge "1" ]; then 
high_min_coarse=`sort -n -k1 Output/period_${period}/2nd_${sta}_coarse | head -n 1 | gawk '{print$2}'`
high_max_coarse=`sort -n -k1 Output/period_${period}/2nd_${sta}_coarse | tail -n 1 | gawk '{print$2}'`
high_source_slow_min=`sort -n -k1 Output/period_${period}/2nd_${sta}_coarse | head -n 1 | gawk '{print$1-0.1}'`
high_source_slow_max=`sort -n -k1 Output/period_${period}/2nd_${sta}_coarse | tail -n 1 | gawk '{print$2+0.1}'`
high_rec_slow_min=`sort -n -k2 Output/period_${period}/2nd_${sta}_coarse | head -n 1 | gawk '{print$2-0.1}'`
high_rec_slow_max=`sort -n -k2 Output/period_${period}/2nd_${sta}_coarse | tail -n 1 | gawk '{print$2+0.1}'`

R1=120.0/121.7/23.5/24.0
R2=0.0/150.0/-1.2/1.2
R3_fund=0.0/150.0/$fund_rec_slow_min/$fund_rec_slow_max
R3_high=0.0/150.0/$high_rec_slow_min/$high_rec_slow_max
R4=0.0/150.0/$dist_min/$dist_max
R5=$high_source_slow_min/$high_source_slow_max/$high_rec_slow_min/$high_rec_slow_max
R6=$fund_source_slow_min/$fund_source_slow_max/$fund_rec_slow_min/$fund_rec_slow_max
J1=M5.5i
J2=X3i/3.5i
J3=X2.6i/1.5i
J4=X2.6i/-2.5i
J5=X2.5i/2.5i
topo_fn=~/Research3/cnliu/topo.grd
topo_fn2=~/Research3/cnliu/taiwan40m.grd
PS=${sta}_per${period}'.ps'
src_sta=`echo $sta | gawk -F_ '{print$1}'`
rec_sta=`echo $sta | gawk -F_ '{print$2}'`
gawk '{print $2,$3}' db_period_${period}/${src_sta}/${rec_sta}/receiver_list > array.file
gawk '{print $2,$3}' station_sort.lst > all.file
gawk '$4==0.0{print $2,$3}' db_period_${period}/${src_sta}/${rec_sta}/receiver_list > center.file
#gawk '$1==sta{print $2,$3}' sta=${src_sta} station_sort.lst > eq.file
gawk '{print $2,$3}' db_period_${period}/${src_sta}/source_list > eq.file
gmt makecpt -Cgray -T-6000/4000/100 -Z > topo.cpt
gawk '{print $1,$2}' Output/period_${period}/1st_${sta}_env > env_f.file
gawk '{print $1,$2}' Output/period_${period}/2nd_${sta}_env > env_h.file

################## Base Map & Station & Array
gmt psbasemap -R$R1 -J$J1 -Bxa0.4f0.2+l"Lon" -Bya0.2f0.1+l"Lat" -BWSen -Gwhite -K > $PS
gmt psxy all.file -R -J -St0.2 -W0.5p,black -O -K>>$PS
gmt psxy array.file -R -J -St0.2 -W0.1p,red -Gred -O -K>>$PS
gmt psxy center.file -R -J -Ss0.3 -W0.1p,black -Gred -O -K>>$PS
gmt psxy eq.file -R -J -Sa0.2 -W0.5p,black -Gyellow -O -K>>$PS
gmt pscoast -R -J -W0.5,black -Df -Swhite -P -K -O >>$PS

################## Stacking Waveform
gmt psbasemap -R$R2 -J$J3 -Bxa20f5+l"Time" -Bya0.5f0.1+l"Amp" -BWSen -Gwhite -Y2.5i -K -O>>$PS
gmt psxy Output/period_${period}/2nd_${sta}_wave -R -J -W0.5p,black -O -K>>$PS
gmt psxy env_h.file -W0.8p,red -R -J -K -O >>$PS
gmt psbasemap -R$R2 -J$J3 -Bxa20f5+l"Time" -Bya0.5f0.1+l"Amp" -BwSen -X3i -Gwhite -K -O>>$PS
gmt psxy Output/period_${period}/1st_${sta}_wave -R -J -W0.5p,black -O -K>>$PS
gmt psxy env_f.file -W0.8p,blue -R -J -K -O >>$PS
################## High-Mode Shift
gmt psbasemap -R$R4 -J$J4 -Bxa20f5+l"Time" -Bya5f1+l"Dist (km)" -BWsen -Gwhite -X-3i -Y2i -K -O>>$PS
gmt pssac Output/period_${period}/2nd_${sta}_shift -R -J -Ekt -M0.9 -Fr -W0.5p,black -O -K>>$PS
################## Fund-Mode Shift
gmt psbasemap -R$R4 -J$J4 -Bxa20f5+l"Time" -Bya5f1+l"Dist (km)" -Bwsen -Gwhite -X3i -K -O>>$PS
gmt pssac Output/period_${period}/1st_${sta}_shift -R -J -Ekt -M0.9 -Fr -W0.5p,black -O -K>>$PS
################## H-Mode Beanforming
max_env_1=`sort -n -k3 Output/period_${period}/2nd_${sta}_slow | tail -n 1 | gawk '{print$3}'`
gawk '{print $1,$2,$3/max_env}' max_env=${max_env_1} Output/period_${period}/2nd_${sta}_slow > slow.file
gmt psbasemap -R$R3_high -J$J2 -Bxa20f5+l"Time" -Bya0.5f0.1+l"Slowness (s/km)" -BWSen -P -K -O -Gwhite -X3.5i -Y-1.0i >> $PS
gmt makecpt -Cjet -T0/1/0.001 -Z > slow.cpt
gmt surface slow.file -R -I0.2/0.01 -Gslow.grd
gmt grdimage slow.grd -Cslow.cpt -R -J -K -O >> $PS
################## Fund-Mode Beanforming
max_env_2=`sort -n -k3 Output/period_${period}/1st_${sta}_slow | tail -n 1 | gawk '{print$3}'`
gawk '{print $1,$2,$3/max_env}' max_env=${max_env_2} Output/period_${period}/1st_${sta}_slow > slow.file
gmt psbasemap -R$R3_fund -J$J2 -Bxa20f5+l"Time" -Bya0.5f0.1+l"Slowness (s/km)" -BwSen -P -K -O -Gwhite -X3.5i >> $PS
gmt makecpt -Cjet -T0/1/0.001 -Z > slow.cpt
gmt surface slow.file -R -I0.2/0.01 -Gslow.grd
gmt grdimage slow.grd -Cslow.cpt -R -J -K -O >> $PS

################## Fund Source-Receiver-Env
max_env=`sort -n -k3 Output/period_${period}/1st_${sta}_coarse | tail -n 1 | gawk '{print$3}'`
gawk '{print $1,$2,$3/max_env}' max_env=${max_env} Output/period_${period}/1st_${sta}_coarse > slow.file
gmt psbasemap -R$R6 -J$J5 -Bxa0.5f0.1+l"Src" -Bya0.5f0.1+l"Rec" -BWSen -P -K -O -Gwhite -Y-3.5i >> $PS
gmt makecpt -Cjet -T0/1/0.002 -Z > slow.cpt
gmt surface slow.file -R -I0.02/0.02 -Gslow.grd
gmt grdimage slow.grd -Cslow.cpt -R -J -K -O >> $PS
################## High Source-Receiver-Env
max_env=`sort -n -k3 Output/period_${period}/2nd_${sta}_coarse | tail -n 1 | gawk '{print$3}'`
gawk '{print $1,$2,$3/max_env}' max_env=${max_env} Output/period_${period}/2nd_${sta}_coarse > slow.file
gmt psbasemap -R$R5 -J$J5 -Bxa0.5f0.1+l"Src" -Bya0.5f0.1+l"Rec" -BWSen -P -K -O -Gwhite -X-3.5i >> $PS
gmt makecpt -Cjet -T0/1/0.002 -Z > slow.cpt
gmt surface slow.file -R -I0.02/0.02 -Gslow.grd
gmt grdimage slow.grd -Cslow.cpt -R -J -K -O >> $PS

else

R1=120.0/121.7/23.5/24.0
R2=0.0/150.0/-1.2/1.2
R3=0.0/150.0/$fund_rec_slow_min/$fund_rec_slow_max
R4=0.0/150.0/$dist_min/$dist_max
R6=$fund_source_slow_min/$fund_source_slow_max/$fund_rec_slow_min/$fund_rec_slow_max
J1=M5.5i
J2=X3i/3.5i
J3=X2.6i/1.5i
J4=X2.6i/-2.5i
J5=X2.5i/2.5i
topo_fn=~/Research3/cnliu/topo.grd
topo_fn2=~/Research3/cnliu/taiwan40m.grd
PS=${sta}_per${period}'.ps'
src_sta=`echo $sta | gawk -F_ '{print$1}'`
rec_sta=`echo $sta | gawk -F_ '{print$2}'`
gawk '{print $2,$3}' db_period_${period}/${src_sta}/${rec_sta}/receiver_list > array.file
gawk '{print $2,$3}' station_sort.lst > all.file
gawk '$4==0.0{print $2,$3}' db_period_${period}/${src_sta}/${rec_sta}/receiver_list > center.file
gawk '{print $2,$3}' db_period_${period}/${src_sta}/source_list > eq.file
makecpt -Cgray -T-6000/4000/100 -Z > topo.cpt
makecpt -Cjet -T0/1/0.01 -Z > slow.cpt
gawk '{print $1,$2}' Output/period_${period}/1st_${sta}_env > env_f.file

################## Base Map & Station & Array
gmt psbasemap -R$R1 -J$J1 -BWSen -Gwhite -K >$PS
gmt psxy all.file -R -J -St0.2 -W0.5p,black -O -K>>$PS
gmt psxy array.file -R -J -St0.2 -W0.1p,red -Gred -O -K>>$PS
gmt psxy center.file -R -J -Ss0.3 -W0.1p,black -Gred -O -K>>$PS
gmt psxy eq.file -R -J -Sa0.2 -W0.5p,black -Gyellow -O -K>>$PS
gmt pscoast -R$R1 -J$J1 -W0.5,black -Df -Swhite -P -K -O>>$PS

################## Stacking Waveform
gmt psbasemap -R$R2 -J$J3 -B50/20WSen -Gwhite -X3i -Y2.5i -K -O>>$PS
gmt psxy Output/period_${period}/1st_${sta}_wave -R -J -W0.5p,black -O -K>>$PS
gmt psxy env_f.file -W0.8p,blue -R -J -K -O >>$PS

################## High-Mode Shift
gmt psbasemap -R$R4 -J$J4 -BWSen -Gwhite -X-3i -Y2i -K -O>>$PS
gmt pssac Output/period_${period}/1st_${sta}_shift -R -J -Ekt -M0.9 -Fr -W0.5p,black -O -K>>$PS

################## Fund-Mode Beanforming
#python src/sort_env.py Output/period_${period}/1st_${sta}_slow
max_env=`sort -n -k3 Output/period_${period}/1st_${sta}_slow | tail -n 1 | gawk '{print$3}'`
gawk '{print $1,$2,$3/max_env}' max_env=${max_env} Output/period_${period}/1st_${sta}_slow > slow.file
gawk '{print $1,$2,$3}' Output/period_${period}/1st_${sta}_slow_sort.txt > slow.file
gmt psbasemap -R$R3 -J$J2 -B50/0.5WSen -P -K -O -Gwhite -X3.5i -Y-1i >> $PS
gmt makecpt -Cjet -T0/1/0.002 -Z > slow.cpt
gmt surface slow.file -R -I0.2/0.01 -Gslow.grd
gmt grdimage slow.grd -Cslow.cpt -R -J -K -O >> $PS
#exit
################## Fund Source-Receiver-Env
max_env=`sort -n -k3 Output/period_${period}/1st_${sta}_coarse | tail -n 1 | gawk '{print$3}'`
gawk '{print $1,$2,$3/max_env}' max_env=${max_env} Output/period_${period}/1st_${sta}_coarse > slow.file
gmt psbasemap -R$R6 -J$J5 -B0.5/0.5WSen -P -K -O -Gwhite -Y-3.5i >> $PS
gmt makecpt -Cjet -T0/1/0.002 -Z > slow.cpt
gmt surface slow.file -R -I0.01/0.01 -Gslow.grd
gmt grdimage slow.grd -Cslow.cpt -R -J -K -O >> $PS
fi
rm *cpt *file *grd *conf *history
convert -density 300 -rotate 90 -trim $PS Figure/$PS'.jpg'
rm $PS
