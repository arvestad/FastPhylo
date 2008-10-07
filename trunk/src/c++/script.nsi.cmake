!include "MUI2.nsh"
Name "fastphylo-${PACKAGE_VERSION}"
outfile "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe"
installDir $PROGRAMFILES\fastphylo

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license-win32-statically-built"
; !insertmacro MUI_PAGE_COMPONENTS 
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
 
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
  
!insertmacro MUI_LANGUAGE "English"

Section "Dummy Section" SecDummy
  SetOutPath "$INSTDIR"
  file fastdist.exe
  file fnj.exe
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\fastdist.exe"
  Delete "$INSTDIR\fnj.exe"
  RMDir "$INSTDIR"
SectionEnd