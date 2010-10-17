; Include files.
!include MUI.nsh
!include AddToPath.nsh
!include FileAssoc.nsh

; Name and output for installation file.
Name "IRTK"
OutFile "IRTKInstaller.exe"

; Default installation folder.
InstallDir "$PROGRAMFILES\IRTK"
  
; Get installation folder from registry if available.
InstallDirRegKey HKLM "Software\IRTK" "InstallationDirectory"

; Interface settings.
!define MUI_ABORTWARNING

; Install pages.
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "SoftwareLicense.txt"
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

; Uninstall pages.
!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

; Languages.
!insertmacro MUI_LANGUAGE "English"

; Main installation section.
Section "Install IRTK" SecInstall
  ; Specfiy read-only section.
  SectionIn RO

  ; Set output path to the installation directory.
  SetOutPath $INSTDIR

  ; Put files there.
  File "..\Windows\bin\release\rview.exe"
  File "..\Windows\bin\release\areg.exe"
  File "..\Windows\bin\release\rreg.exe"
  File "..\Windows\bin\release\nreg.exe"
  File "..\Windows\bin\release\pareg.exe"
  File "..\Windows\bin\release\phreg.exe"
  File "..\Windows\bin\release\pnreg.exe"
  File "..\Windows\bin\release\prreg.exe"
  File "..\Windows\bin\release\sareg.exe"
  File "..\Windows\bin\release\shreg.exe"
  File "..\Windows\bin\release\snreg.exe"
  File "..\Windows\bin\release\srreg.exe"
  File "..\Windows\bin\release\motiontrack.exe"
  File "..\Windows\bin\release\transformation.exe"
  File "..\Windows\bin\release\stransformation.exe"
  File "..\Windows\bin\release\ptransformation.exe"
  File "..\Windows\bin\release\jacobian.exe"
  File "..\Windows\bin\release\atlas.exe"
  File "..\Windows\bin\release\dmap.exe"
  File "..\Windows\bin\release\dof2flirt.exe"
  File "..\Windows\bin\release\dof2image.exe"
  File "..\Windows\bin\release\dof2mat.exe"
  File "..\Windows\bin\release\dofinvert.exe"
  File "..\Windows\bin\release\evaluation.exe"
  File "..\Windows\bin\release\flirt2dof.exe"
  File "..\Windows\bin\release\info.exe"
  File "..\Windows\bin\release\convert.exe"
  File "..\Windows\bin\release\threshold.exe"
  File "..\Windows\bin\release\padding.exe"
  File "..\Windows\bin\release\blur.exe"
  File "..\Windows\bin\release\dilation.exe"
  File "..\Windows\bin\release\dmap.exe"
  File "..\Windows\bin\release\makesequence.exe"
  File "..\Windows\bin\release\erosion.exe"
  File "..\Windows\bin\release\closing.exe"
  File "..\Windows\bin\release\opening.exe"
  File "..\Windows\bin\release\reflect.exe"
  File "..\Windows\bin\release\region.exe"
  File "..\Windows\bin\release\resample.exe"
  File "..\Windows\bin\release\rescale.exe"
  File "..\README"
  File "..\COPYRIGHT"

  File "..\icons\IRTK.ico"

  ;File "C:\Program Files\Microsoft Visual Studio 8\VC\redist\x86\Microsoft.VC80.CRT\msvcp80.dll"
  ;File "C:\Program Files\Microsoft Visual Studio 8\VC\redist\x86\Microsoft.VC80.CRT\msvcr80.dll"
 
  ; Write the installation path into the registry.
  WriteRegStr HKLM SOFTWARE\IRTK "InstallationDirectory" "$INSTDIR"

  ; Write the uninstall keys for Windows.
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\IRTK" "DisplayName" "IRTK"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\IRTK" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\IRTK" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\IRTK" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
SectionEnd

; Optional section for creating start menu shortcuts.
Section "Start menu shortcuts" SecStartMenuShortcuts
  CreateDirectory "$SMPROGRAMS\IRTK"
  CreateShortCut "$SMPROGRAMS\IRTK\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\IRTK\rview.lnk" "$INSTDIR\rview.exe" "" "$INSTDIR\rview.exe" "0"
SectionEnd

; Optional section for adding directories to the PATH environment variable.
Section "Environment Variables" SecEnvironmentVariables
  Push $INSTDIR
  Call AddToPath
SectionEnd

; Optional section for associating .hdr, .gipl, and .nii files
; with rview.
Section "File Associations" SecFileAssociations
  !insertmacro APP_ASSOCIATE "hdr" "rview.hdrFile" "$INSTDIR\rview.exe" "$INSTDIR\IRTK.ico" "Open with rview" "$INSTDIR\rview.exe $\"%1$\""
  !insertmacro APP_ASSOCIATE "gipl" "rview.giplFile" "$INSTDIR\rview.exe" "$INSTDIR\IRTK.ico" "Open with rview" "$INSTDIR\rview.exe $\"%1$\""
  !insertmacro APP_ASSOCIATE "nii" "rview.niiFile" "$INSTDIR\rview.exe" "$INSTDIR\IRTK.ico" "Open with rview" "$INSTDIR\rview.exe $\"%1$\""
SectionEnd

; Descriptions.
LangString DESC_SecInstall ${LANG_ENGLISH} "Install Image Registration Toolkit (IRTK) software."
LangString DESC_SecStartMenuShortcuts ${LANG_ENGLISH} "Add shortcuts to frequently used programs to the Windows Start Menu."
LangString DESC_SecEnvironmentVariables ${LANG_ENGLISH} "Add the installation directory to the PATH environment variable so that IRTK programs can be executed from the command prompt without qualifying them with the full path name."
LangString DESC_SecFileAssociations ${LANG_ENGLISH} "Associate .hdr and .gipl files with rview."

; Assign language strings to sections.
!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
  !insertmacro MUI_DESCRIPTION_TEXT ${SecInstall} $(DESC_SecInstall)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecStartMenuShortcuts} $(DESC_SecStartMenuShortcuts)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecEnvironmentVariables} $(DESC_SecEnvironmentVariables)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecFileAssociations} $(DESC_SecFileAssociations)
!insertmacro MUI_FUNCTION_DESCRIPTION_END

; Uninstaller section.
Section "Uninstall"
  ; Remove the registry keys.
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\IRTK"
  DeleteRegKey HKLM "Software\IRTK"

  ; Remove files, uninstaller, and icon
  Delete $INSTDIR\*.exe
  Delete $INSTDIR\IRTK.ico

  ; Remove shortcuts.
  Delete "$SMPROGRAMS\IRTK\*.*"

  ; Remove directories used.
  RMDir "$SMPROGRAMS\IRTK"
  RMDir "$INSTDIR"

  ; Remove installation directory from path.
  Push $INSTDIR
  Call un.RemoveFromPath

  ; Remove file associations for rview.
  !insertmacro APP_UNASSOCIATE "hdr" "rview.hdrFile"
  !insertmacro APP_UNASSOCIATE "gipl" "rview.giplFile"
  !insertmacro APP_UNASSOCIATE "nii" "rview.niiFile"
SectionEnd
