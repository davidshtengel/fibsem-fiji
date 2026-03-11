@echo off	
rem Set FIJI_PLUGINS to the path of Fiji\plugins on your machine 
set "FIJI_PLUGINS=C:\Users\User\Fiji\plugins"

echo Building ...
rem Use Maven to remove files from previous build, then compile and package into .jar
call mvn clean package -q 
if errorlevel 1 (
    echo ERROR: Build failed.
    exit /b 1
)

echo Removing old plugin JARs ...
rem Use pattern-matching to removing any existing/older versions of the plugin
del "%FIJI_PLUGINS%\fib_sem-*.jar" 2>nul

echo Installing ...
rem Use findstr /r and regex to copy just the fib_sem-#.#.#.jar into the Fiji\plugins directory specified above
for %%f in (target\fib_sem-*.jar) do (
    echo %%~nxf | findstr /r "^fib_sem-[0-9][0-9]*[.][0-9][0-9]*[.][0-9][0-9]*[.]jar$" >nul && copy /y "%%f" "%FIJI_PLUGINS%\" >nul
)

echo Done! Restart Fiji.
pause