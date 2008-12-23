#!/bin/csh

set targetExe = "rreg areg nreg pareg pnreg prreg sareg snreg srreg motiontrack transformation stransformation ptransformation jacobian atlas dmap dof2flirt dof2image dof2mat dofinvert evaluation flirt2dof info convert threshold binarize mcubes padding blur dilation dmap erosion closing makesequence opening reflect region resample rescale rview"

set targetDir = "itk"

mkdir -p $targetDir

foreach f ( $targetExe )
   
   echo "Copying $f to directory $targetDir"
   cp ${ITK_BINARY_DIR}/bin/$f $targetDir

end

cp ${ITK_SOURCE_DIR}/README $targetDir
cp ${ITK_SOURCE_DIR}/COPYRIGHT $targetDir

if ( `uname` == "Linux" ) then

   tar czvf itk-linux-32.tar.gz $targetDir

else

   hdiutil create -srcfolder $targetDir itk-mac-32.dmg

endif
