@echo off
:: Batch script to replace "Physiolibrary." with "Bodylight."
:: while preserving UTF-8 encoding and ignoring "/Physiolibrary."

echo Starting replacement of "Physiolibrary." with "Bodylight." in .mo files...

:: Loop through all .mo files in the current directory and subdirectories
for /r %%f in (*.mo) do (
    echo Processing file: %%f
    powershell -Command "(Get-Content -Raw -Encoding UTF8 '%%f') -replace '(?<!/)Physiolibrary\.', 'Bodylight.' | Set-Content -NoNewline -Encoding UTF8 '%%f'"
)

echo Replacement completed!
pause
