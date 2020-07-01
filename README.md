# raytracer

HOW TO RUN:
1. Turn on command line on windows 10 (type "cmd" in start menu) and go to the directory you extracted this zip file (type "cd pathname" on command line).
2. Go to "raytracer executable" folder ("cd raytracer executable"), then type like "raytracer.exe testAmbient.txt" then press enter.
3. It will generate ppm file per each command you run with txt file in same directory.

EXTRA INFO:
Coded the program on Visual Studio 2017 Community version on Windows 10.
If you just run the program by typing "raytracer.exe testCase1.txt" in command line without Visual Studio 2017 installed on your machine, normally it will give you errors if you do not have Visual Studio 2017 installed on your machine.
So you will need to have 2 DLL files to run this program, 32bit version of both "msvcp140.dll" and "vcruntime140.dll". Fortunately, I have put them all in "raytracer executable" directory along with "raytracer.exe" file.
Thus, you should not have a problem to run this program if you are on same WINDOWS 10. Otherwise, please install Visual Studio and open sln file in "raytracer source" directory, compile in release mode and run manually.

![alt text](https://i.imgur.com/iHJNvRp.jpg)

![alt text](https://i.imgur.com/gmAxJJf.jpg)
