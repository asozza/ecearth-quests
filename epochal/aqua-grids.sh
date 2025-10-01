
#!/bin/bash
DIR=/lus/h2resw01/hpcperm/ccpd/AQUAGRID

# cleanup
module load cdo/2.4.0

rm -rf $DIR
mkdir -p $DIR

#ORCA2
cdo="cdo -f nc1 -b 32"

#area
$cdo -expr,areat=e1t*e2t /lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain/ORCA2/domain_cfg.nc $DIR/cellarea_ORCA2.nc

# ORCA2
exp=TN10
user=ccpd
mkdir -p $DIR/ORCA2
$cdo gtc,0 -setname,mask -selname,sos -seltimestep,1 /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/ORCA2/ORCA2_mesh_sfc_grid_T.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,thetao /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/ORCA2/ORCA2_mesh_3d_grid_T.nc
#rm -f $DIR/tmp.nc

#PALEORCA

#area
$cdo -expr,areat=e1t*e2t /lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain/PALEORCA2/domain_cfg.nc $DIR/cellarea_PALEORCA2.nc

# PALEORCA-PD
exp=LC05
user=ecme3497
mkdir -p $DIR/PALEORCA2
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,sos /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/PALEORCA2_mesh_sfc_grid_T-pd.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,so /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/PALEORCA2_mesh_3d_grid_T-pd.nc



# PALEORCA-EO
exp=EP01
user=ecme3497
mkdir -p $DIR/PALEORCA2
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,sos /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/PALEORCA2_mesh_sfc_grid_T-eo.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,thetao /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/PALEORCA2_mesh_3d_grid_T-eo.nc