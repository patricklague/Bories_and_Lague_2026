aafile=('SCW' 'SCF' 'SCY')

echo "" > out.out

for aa in "${aafile[@]}"
do
  echo $aa
  python ring_orientation_analysis.py $aa
  python freq_angle_analysis.py $aa
  bash ./extract_top_angles.sh $aa >> out.out
done
mv orientation_* ../../data/aromatics_orientation/raw_data/
mv freq_angle_* ../../data/aromatics_orientation/
mv out.out ../../data/aromatics_orientation/top_angles.out
