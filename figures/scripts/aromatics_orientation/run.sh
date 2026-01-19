aafile=('SCW' 'SCF' 'SCY')

echo "" > out_mono.out

for aa in "${aafile[@]}"
do
  echo $aa
  python ring_orientation_analysis.py $aa
  python freq_angle_analysis.py $aa
  bash ./extract_top_angles.sh $aa >> out_mono.out
done
mv orientation_mono_* ../../data/aromatics_orientation/raw_data/
mv freq_angle_mono_* ../../data/aromatics_orientation/
mv out_mono.out ../../data/aromatics_orientation/top_angles_mono.out
