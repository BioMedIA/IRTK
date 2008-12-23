; =========================================================================
; 
;  Library   : ITK Installer
;  Module    : $RCSfile: itkInstaller.nsi,v $
;  Authors   : Raghavendra Chandrashekara
;  Copyright : Imperial College, Department of Computing
;              Visual Information Processing (VIP), 2006
;  Purpose   : Installer script for ITK.
;  Date      : $Date: 2008-12-05 16:14:00 $
;  Version   : $Revision: 1.12 $
;  Changes   : $Locker:  $
;              $Log: itkInstaller.nsi,v $
;              Revision 1.12  2008-12-05 16:14:00  dr
;              *** empty log message ***
;
;              Revision 1.11  2008-12-01 10:02:36  dr
;              *** empty log message ***
;
;              Revision 1.10  2008-10-31 12:25:41  dr
;              Added new executables
;
;              Revision 1.9  2008-05-25 10:26:05  dr
;              Fixed another spelling mistake
;
;              Revision 1.8  2008-05-25 10:24:06  dr
;              Fixed spelling mistake
;
;              Revision 1.7  2008-05-24 17:39:05  dr
;              *** empty log message ***
;
;              Revision 1.6  2007-12-17 21:01:41  dr
;              *** empty log message ***
;
;              Revision 1.5  2006-08-23 21:26:36  dr
;              Added support for NIFTI
;
;              Revision 1.4  2006/08/15 21:40:18  dr
;              *** empty log message ***
;
;              Revision 1.3  2006/03/28 17:23:39  rc3
;              Added license page.
;
;              Revision 1.2  2006/03/27 17:27:27  rc3
;              Updated script to use the "Modern User Interface"
;              (MUI) macros.
;
;              Revision 1.1  2006/03/27 16:14:06  rc3
;              Added scripts for installing ITK software on Windows.
;
; =========================================================================

; Include files.
!include MUI.nsh
!include AddToPath.nsh
!include FileAssoc.nsh

; Name and output for installation file.
Name "ITK"
OutFile "ITKInstaller.exe"

; Default installation folder.
InstallDir "$PROGRAMFILES\ITK"
  
; Get installation folder from registry if available.
InstallDirRegKey HKLM "Software\ITK" "InstallationDirectory"

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
Section "Install ITK" SecInstall
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

  File "..\icons\ITK.ico"

  File "C:\Program Files\Microsoft Visual Studio 8\VC\redist\x86\Microsoft.VC80.CRT\msvcp80.dll"
  File "C:\Program Files\Microsoft Visual Studio 8\VC\redist\x86\Microsoft.VC80.CRT\msvcr80.dll"
 
  ; Write the installation path into the registry.
  WriteRegStr HKLM SOFTWARE\ITK "InstallationDirectory" "$INSTDIR"

  ; Write the uninstall keys for Windows.
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\ITK" "DisplayName" "ITK"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\ITK" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\ITK" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\ITK" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
SectionEnd

; Optional section for creating start menu shortcuts.
Section "Start menu shortcuts" SecStartMenuShortcuts
  CreateDirectory "$SMPROGRAMS\ITK"
  CreateShortCut "$SMPROGRAMS\ITK\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\ITK\rview.lnk" "$INSTDIR\rview.exe" "" "$INSTDIR\rview.exe" "0"
SectionEnd

; Optional section for adding directories to the PATH environment variable.
Section "Environment Variables" SecEnvironmentVariables
  Push $INSTDIR
  Call AddToPath
SectionEnd

; Optional section for associating .hdr, .gipl, and .nii files
; with rview.
Section "File Associations" SecFileAssociations
  !insertmacro APP_ASSOCIATE "hdr" "rview.hdrFile" "$INSTDIR\rview.exe" "$INSTDIR\ITK.ico" "Open with rview" "$INSTDIR\rview.exe $\"%1$\""
  !insertmacro APP_ASSOCIATE "gipl" "rview.giplFile" "$INSTDIR\rview.exe" "$INSTDIR\ITK.ico" "Open with rview" "$INSTDIR\rview.exe $\"%1$\""
  !insertmacro APP_ASSOCIATE "nii" "rview.niiFile" "$INSTDIR\rview.exe" "$INSTDIR\ITK.ico" "Open with rview" "$INSTDIR\rview.exe $\"%1$\""
SectionEnd

; Descriptions.
LangString DESC_SecInstall ${LANG_ENGLISH} "Install Image Registration Toolkit (ITK) software."
LangString DESC_SecStartMenuShortcuts ${LANG_ENGLISH} "Add shortcuts to frequently used programs to the Windows Start Menu."
LangString DESC_SecEnvironmentVariables ${LANG_ENGLISH} "Add the installation directory to the PATH environment variable so that ITK programs can be executed from the command prompt without qualifying them with the full path name."
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
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\ITK"
  DeleteRegKey HKLM "Software\ITK"

  ; Remove files, uninstaller, and icon
  Delete $INSTDIR\*.exe
  Delete $INSTDIR\ITK.ico

  ; Remove shortcuts.
  Delete "$SMPROGRAMS\ITK\*.*"

  ; Remove directories used.
  RMDir "$SMPROGRAMS\ITK"
  RMDir "$INSTDIR"

  ; Remove installation directory from path.
  Push $INSTDIR
  Call un.RemoveFromPath

  ; Remove file associations for rview.
  !insertmacro APP_UNASSOCIATE "hdr" "rview.hdrFile"
  !insertmacro APP_UNASSOCIATE "gipl" "rview.giplFile"
  !insertmacro APP_UNASSOCIATE "nii" "rview.niiFile"
SectionEnd
