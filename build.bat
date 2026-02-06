@echo off
setlocal

:: Set up MSVC environment
call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64

if not exist build mkdir build
cd build

cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 (
    echo CMake configuration failed.
    exit /b 1
)

cmake --build .
if errorlevel 1 (
    echo Build failed.
    exit /b 1
)

echo.
echo Build successful! Executable: build\MeshShaderFont.exe
