#!/bin/csh

set targetExe = "ems combineLabels labelStats rreg areg nreg pareg pnreg prreg sareg snreg srreg motiontrack transformation stransformation ptransformation jacobian atlas dmap dof2flirt dof2image dof2mat dofinvert evaluation flirt2dof info image2dof convert threshold binarize mcubes padding blur dilation dmap erosion closing makesequence opening reflect region resample rescale rview"

set targetDir = "irtk"

mkdir -p $targetDir

foreach f ( $targetExe )
   
   echo "Copying $f to directory $targetDir"
   cp ${IRTK_BINARY_DIR}/bin/$f $targetDir

end

cp ${IRTK_SOURCE_DIR}/README $targetDir
cp ${IRTK_SOURCE_DIR}/COPYRIGHT $targetDir

if ( `uname` == "Linux" ) then

   tar czvf irtk-linux-64.tar.gz $targetDir

else

   hdiutil create -srcfolder $targetDir irtk-mac-64.dmg

endif
