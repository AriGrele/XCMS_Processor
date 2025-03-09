if not defined in_subprocess (cmd /k set in_subprocess=y ^& %0 %*) & exit )
@echo off
set /P mes=Message? 
git add .
git commit -m %mes%
git push origin HEAD:master