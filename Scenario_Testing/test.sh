
cd .. 

echo "Testing Homogeneous model"
pandemia Scenarios/Homogeneous/homogeneous_config.yaml

echo "Testing Heterogeneous model"
pandemia Scenarios/Heterogeneous/heterogeneous_config.yaml

echo "Testing saving world file"
pandemia Scenarios/Heterogeneous/heterogeneous_config.yaml Scenarios/Heterogeneous/heterogeneous_world.wld

echo "Testing saved world with mixing"
pandemia Scenarios/Heterogeneous/heterogeneous_config.yaml Scenarios/Heterogeneous/heterogeneous_world.wld

cd Scenario_Testing