
#!/bin/bash
DIR=/ec/res4/hpcperm/itmn/AQUAGRID

# cleanup
module load cdo/2.4.0
module load nco/5.3.2

rm -rf $DIR
mkdir -p $DIR

cdo="cdo -f nc1 -b 32"

#----ORCA2------#
exp=TN10
user=ccpd
mkdir -p $DIR/ORCA2
$cdo gtc,0 -setname,mask -selname,sos -seltimestep,1 /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/ORCA2/ORCA2_mesh_sfc_grid_T.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,thetao /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/ORCA2/ORCA2_mesh_3d_grid_T.nc
for file in $DIR/ORCA2/ORCA2_mesh_sfc_grid_T.nc $DIR/ORCA2/ORCA2_mesh_3d_grid_T.nc ; do
  ncrename -d x_grid_T,x $file
  ncrename -d y_grid_T,y $file
done

#area
$cdo -expr,area=e1t*e2t /lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain/ORCA2/domain_cfg.nc $DIR/ORCA2/cellarea_ORCA2_T.nc

#----PALEORCA-----#

# PALEORCA-PD
exp=LC05
user=ecme3497
grid=PALEORCA2
window=pd
mkdir -p $DIR/$grid
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,sos /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/${grid}_mesh_sfc_grid_T-${window}.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,so /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/${grid}_mesh_3d_grid_T-${window}.nc

# PALEORCA-EO
exp=EP01
user=ecme3497
window=eo
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,sos /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/${grid}_mesh_sfc_grid_T-${window}.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,thetao /lus/h2resw01/scratch/${user}/ece4/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/PALEORCA2/${grid}_mesh_3d_grid_T-${window}.nc

for window in pd eo ; do
  for file in sfc 3d ; do
    ncrename -d x_grid_T,x $DIR/PALEORCA2/${grid}_mesh_${file}_grid_T-${window}.nc
    ncrename -d y_grid_T,y $DIR/PALEORCA2/${grid}_mesh_${file}_grid_T-${window}.nc
  done
done

#area
$cdo -expr,area=e1t*e2t /lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain/PALEORCA2/domain_cfg.nc $DIR/PALEORCA2/cellarea_PALEORCA2_T.nc


#-----eORCA1-----#
exp=PR04
user=smw
mkdir -p $DIR/eORCA1
$cdo gtc,0 -setname,mask -selname,sos -seltimestep,1 /perm/${user}/ece-4-exps/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/eORCA1/eORCA1_mesh_sfc_grid_T.nc
$cdo gtc,0 -setname,mask -seltimestep,1 -selname,thetao /perm/${user}/ece-4-exps/${exp}/output/nemo/${exp}_oce_1m_T_1990-1990.nc $DIR/eORCA1/eORCA1_mesh_3d_grid_T.nc
for file in $DIR/eORCA1/eORCA1_mesh_sfc_grid_T.nc $DIR/eORCA1/eORCA1_mesh_3d_grid_T.nc ; do
  ncrename -d x_grid_T,x $file
  ncrename -d y_grid_T,y $file
done

#area
$cdo -expr,area=e1t*e2t /lus/h2resw01/hpcperm/ccpd/ECE4-DATA/nemo/domain/eORCA1/domain_cfg.nc $DIR/eORCA1/cellarea_eORCA1_T.nc
