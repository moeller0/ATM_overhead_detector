@ECHO OFF

SETLOCAL

REM this could be smarter...

REM try to run hrping to collect data to investigate

REM whether a link uses ATM cell quantization




REM Is hrping installed? Bail out and complain if not

REM Is hrping installed? Bail out & complain if not

REM Get hrping from: https://www.cfos.de/en/ping/ping.htm
REM hrping required EULA confirmation on initial run

ECHO This script requires administrative priviledge


REM Define variables
REM gstatic.com does seem to respond only up to 64 bytes ICMP payload...

SET TARGET=8.8.8.8

SET /A SWEEPMINSIZE=16
SET /A SWEEP_N_ATM_CELLS=3
SET /A SWEEPMAXSIZE=%SWEEP_N_ATM_CELLS% * 48 + %SWEEPMINSIZE%

SET PINGPERIOD_MS=10
SET PINGSPERSIZE=1000


SET /A N_SIZES=%SWEEPMAXSIZE% - %SWEEPMINSIZE%
ECHO %N_SIZES%
SET /A N_PINGSTOTAL=%N_SIZES% * %PINGSPERSIZE%
ECHO %N_PINGSTOTAL%

SET TECH=xDSL

REM get time and date and remove some unwanted characters
SET TIMESTR=%TIME: =%
SET TIMESTR=%TIMESTR::=%
SET DATESTR=%DATE: =%

SET DATESTR=%DATESTR:/=%

SET LOG_FILE_NAME=ping_sweep_%TECH%_%DATESTR%_%TIMESTR%.txt

SET PING_CMD=hrping -q -l%SWEEPMINSIZE%:%SWEEPMAXSIZE%:1 -s %PINGPERIOD_MS% -W -n %N_PINGSTOTAL% -F %LOG_FILE_NAME% %TARGET%

ECHO %PING_CMD%


REM now run this
%PING_CMD%

ECHO DONE
ENDLOCAL
ECHO ON