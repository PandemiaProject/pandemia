
cd .. 

echo "Testing global model"
pandemia Scenarios/Global/global_config.yaml

echo "Testing global grid model"
pandemia Scenarios/Global_Grid/global_grid_config.yaml

echo "Testing saving world file"
pandemia Scenarios/Global_Grid/global_grid_config.yaml Scenarios/Global_Grid/global_grid_world.wld

echo "Testing saved world with mixing"
pandemia Scenarios/Global_Grid/global_grid_config.yaml Scenarios/Global_Grid/global_grid_world.wld

cd Scenario_Testing