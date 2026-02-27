aafile=('SCW' 'SCF' 'SCY') #'SCF' 'SCY' 'SCYM' 'SCHP' 'SCHE' 'SCHD' 'SCW')

for aa in "${aafile[@]}"
do
  echo $aa
  python ring_orientation_analysis.py $aa
  python freq_angle_analysis.py $aa
  bash ./extract_top_angles.sh $aa >> out_mono.out
done
