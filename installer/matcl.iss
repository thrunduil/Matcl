; Inno setup installer
; http://www.jrsoftware.org/isinfo.php

#include "environment.iss"

[Setup]

AppVersion=1.0.1
#define ver "1.0.1"

AppName=matcl
DefaultDirName={pf}\matcl

; publisher
AppPublisher=matcl
AppPublisherURL=https://bitbucket.org/matcl/matcl/src

; Setup will not show the Select Start Menu Folder wizard page.
DisableProgramGroupPage=yes
DefaultGroupName=matcl

UninstallDisplayIcon={app}\matcl_ico.ico 
UninstallDisplayName=matcl_{#ver}

OutputBaseFilename=matcl_{#ver}_setup

ChangesEnvironment=true

Compression=lzma2
SolidCompression=yes

LicenseFile="..\LICENSE"

; "ArchitecturesAllowed=x64" specifies that Setup cannot run on
; anything but x64.
ArchitecturesAllowed=x64

; "ArchitecturesInstallIn64BitMode=x64" requests that the install be
; done in "64-bit mode" on x64, meaning it should use the native
; 64-bit Program Files directory and the 64-bit view of the registry.
ArchitecturesInstallIn64BitMode=x64

OutputDir="."

SetupIconFile=matcl_ico.ico
WizardImageFile=matcl_bmp.bmp
WizardSmallImageFile=matcl_bmp.bmp

; uninstallers
UninstallFilesDir={app}\{#ver}\uninstall

; hide uninstallers
[Dirs]

Name: {app}\{#ver}\uninstall; Attribs: hidden;

[Files]

; other files
Source: "matcl_ico.ico"; DestDir: "{app}"

; include files
Source: "configs\matcl_config_release.h"; DestDir: "{app}\{#ver}\include\"; DestName: "matcl_config.h";

Source: "..\src\include\matcl-core\*"; DestDir: "{app}\{#ver}\include\matcl-core"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-dynamic\*"; DestDir: "{app}\{#ver}\include\matcl-dynamic"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-file\*"; DestDir: "{app}\{#ver}\include\matcl-file"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-matrep\*"; DestDir: "{app}\{#ver}\include\matcl-matrep"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-mp\*"; DestDir: "{app}\{#ver}\include\matcl-mp"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-mp-obj\*"; DestDir: "{app}\{#ver}\include\matcl-mp-obj"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-scalar\*"; DestDir: "{app}\{#ver}\include\matcl-scalar"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-simd\*"; DestDir: "{app}\{#ver}\include\matcl-simd"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-specfunc\*"; DestDir: "{app}\{#ver}\include\matcl-specfunc"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-internals\*"; DestDir: "{app}\{#ver}\include\matcl-internals"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\src\include\matcl-blas-lapack\*"; DestDir: "{app}\{#ver}\include\matcl-blas-lapack"; Flags: ignoreversion createallsubdirs recursesubdirs

Source: "..\src\include\matcl-linalg\*"; DestDir: "{app}\{#ver}\include\matcl-linalg"; Flags: ignoreversion createallsubdirs recursesubdirs

; lib files 64 bit
Source: "..\x64\release\matcl-matrep-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-matfunc-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-matmult-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-core-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-dynamic-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-scalar-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-file-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-specfunc-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-simd-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-mp-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-mp-obj-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-blas-loader-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-blas-lapack-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";
Source: "..\x64\release\matcl-flapack-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";

Source: "..\x64\release\linalg-x64-Release.lib"; DestDir: "{app}\{#ver}\x64\lib";

; lib files 32 bit
Source: "..\win32\release\matcl-matrep-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-matfunc-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-matmult-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-core-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-dynamic-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-scalar-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-file-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-specfunc-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-simd-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-mp-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-mp-obj-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-blas-loader-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-blas-lapack-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";
Source: "..\win32\release\matcl-flapack-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";

Source: "..\win32\release\linalg-win32-Release.lib"; DestDir: "{app}\{#ver}\win32\lib";

; bin files 64 bit
Source: "..\x64\release\matcl-matmult-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-core-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-scalar-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-dynamic-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-matrep-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-matfunc-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-blas-lapack-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-flapack-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-blas-loader-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-mp-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\libquadmath-0.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\mpir-x64-release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\mpfr-x64-release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\libgcc_s_seh-1.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-file-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-mp-obj-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-sqlite-cpp-x64-release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-specfunc-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-openblas-plugin-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-simd-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-clapack-plugin-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\openblas-x64-release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\libgfortran-4.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\libwinpthread-1.dll"; DestDir: "{app}\{#ver}\x64\bin";

Source: "..\x64\release\blas_config_win64.txt"; DestDir: "{app}\{#ver}\x64\bin";

Source: "..\x64\release\linalg-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl_arpack-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl_metis-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\matcl-blas-ext-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\scotch-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\cholmod-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";
Source: "..\x64\release\superlu-x64-Release.dll"; DestDir: "{app}\{#ver}\x64\bin";

; bin files 32 bit
Source: "..\win32\release\matcl-matmult-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-core-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-scalar-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-dynamic-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-matrep-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-matfunc-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-blas-lapack-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-flapack-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-blas-loader-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-mp-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\mpir-win32-release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\mpfr-win32-release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-file-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-mp-obj-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-sqlite-cpp-win32-release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-specfunc-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-openblas-plugin-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-simd-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-clapack-plugin-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\openblas-win32-release.dll"; DestDir: "{app}\{#ver}\win32\bin";

Source: "..\win32\release\blas_config_win32.txt"; DestDir: "{app}\{#ver}\win32\bin";

Source: "..\win32\release\linalg-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl_arpack-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl_metis-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\matcl-blas-ext-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\scotch-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\cholmod-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";
Source: "..\win32\release\superlu-win32-Release.dll"; DestDir: "{app}\{#ver}\win32\bin";

[Icons]
Name: "{group}\{cm:UninstallProgram,matlc_{#ver}}"; Filename: "{uninstallexe}"; IconFilename: matcl_ico.ico

; shortcut to unistaller
Name: "{app}\matcl_{#ver}_uninstall"; Filename: "{app}\{#ver}\uninstall\unins000.exe"

[Code]
procedure CurStepChanged(CurStep: TSetupStep);
begin
    if CurStep = ssPostInstall 
    then EnvAddPath(ExpandConstant('{app}\{#ver}') +'\x64\bin');
         EnvAddPath(ExpandConstant('{app}\{#ver}') +'\win32\bin');
end;

procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
begin
    if CurUninstallStep = usPostUninstall
    then EnvRemovePath(ExpandConstant('{app}\{#ver}') +'\x64\bin');
         EnvRemovePath(ExpandConstant('{app}\{#ver}') +'\win32\bin');
end;
