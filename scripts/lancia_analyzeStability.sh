for i in {1243..1400};
do 
run=$(printf "%04d" $i)
#sed -i "s%^recoFileNames .*$%recoFileNames run"$run"_ped_t.root%" cfg/analyzeStability.cfg;
#sed -i "s%^rawFileNames .*$%rawFileNames run"$run".rawf%" cfg/analyzeStability.cfg;
#sed -i "s%^outFileName .*$%outFileName /home/cmsdaq/Lab5015Analysis/plots/analyzeStability/run"$run"_ped_stability.root%" cfg/analyzeStability.cfg;
#sed -i "s%^plotDir .*$%plotDir /var/www/html/stability/run"$run"_ped%" cfg/analyzeStability.cfg;
#./bin/analyzeStability.exe cfg/analyzeStability.cfg
sed -i "s%^recoFileNames .*$%recoFileNames run"$run"_t.root%" cfg/analyzeStability.cfg;
sed -i "s%^rawFileNames .*$%rawFileNames run"$run".rawf%" cfg/analyzeStability.cfg;
sed -i "s%^outFileName .*$%outFileName /home/cmsdaq/Lab5015Analysis/plots/analyzeStability/run"$run"_stability.root%" cfg/analyzeStability.cfg;
sed -i "s%^plotDir .*$%plotDir /var/www/html/stability/run"$run"%" cfg/analyzeStability.cfg;
./bin/analyzeStability.exe cfg/analyzeStability.cfg
done;
