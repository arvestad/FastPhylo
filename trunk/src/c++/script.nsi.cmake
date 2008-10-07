!include "MUI2.nsh"
Name "fastphylo"
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

  ;Language strings
;  LangString DESC_SecDummy ${LANG_ENGLISH} "en minitestsektion."

  ;Assign language strings to sections
;  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
;    !insertmacro MUI_DESCRIPTION_TEXT ${SecDummy} $(DESC_SecDummy)
;  !insertmacro MUI_FUNCTION_DESCRIPTION_END

Section "Uninstall"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\fastdist.exe"
  Delete "$INSTDIR\fnj.exe"
  RMDir "$INSTDIR"
SectionEnd