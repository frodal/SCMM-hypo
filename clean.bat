@echo off
if exist HypoExp-xpl.obj del HypoExp-xpl.obj
if exist HypoExp-xplD.obj del HypoExp-xplD.obj
if exist HypoImp-std.obj del HypoImp-std.obj
if exist explicitU.dll del explicitU.dll
if exist explicitU-D.dll del explicitU-D.dll
if exist standardU.dll del standardU.dll
call python Test/Test.py clean