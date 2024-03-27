# SAMPLES=( "10X-Bezerra-CO-022-MB0-20211117-3v3-1hg" "10X-Bezerra-CO-53-MBO-NC-20211201-3v3-1hg" "10X-Bezerra-CO-BLP-72-MBO-20211020-3v3-1hg" "10X-Bezerra-CO21-MBO-20211123-3v3-1hg" "10X-Bezerra-CO36-MBO-20211103-3v3-1hg" "10X-Bezerra-CO55-MBO-20211027-3v3-1hg" )

# SAMPLES=( "10X-Bezerra-CO-68-MBO-20220302-3v3" )

SAMPLES=( "10X-Bezerra-CO72-MBO-BA-p10-20220509" "10X-Bezerra-CO21_MBO-BA-20220621" )

for i in "${SAMPLES[@]}"
do
	echo $i
	bsub -L /bin/bash -W 40:00 -n 4 -R "span[ptile=4]" -M 128000 -e ./${i}.err -o ./${i}.out ./cellranger6.sh ${i}
done
