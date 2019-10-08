#########################################################################
#                                                                       #
#    DGT tracker - software for Robust Deformable and Occluded Object   #
#    Tracking with Dynamic Graph, Version 2.1                           #
#    https://sites.google.com/site/zhaoweicai1989/                      #
#                                                                       #
#    Copyright 2015 Zhaowei Cai (zwcai@ucsd.edu)                        #
#                                                                       #
#########################################################################

/************************************  IMPORTANT!!!!!!!!!!!  ****************************************/

  If you use this software, YOU should cite the following 2 papers IN ANY RESULTING PUBLICATION:
  
  [1] Robust Deformable and Occluded Object Tracking with Dynamic Graph,
			Zhaowei Cai, Longyin Wen, Zhen Lei, Nuno Vasconcelos, Stan Z. Li. IEEE Transaction on Image Processing (TIP), 2014 
  [2] Structured Visual Tracking with Dynamic Graph,
			Zhaowei Cai, Longyin Wen, Jianwei Yang, Zhen Lei, Stan Z. Li. Asian Conference on Computer Vision (ACCV), 2012 

/****************************************************************************************************/
   In order to use this code you should
   a. Add the path and dll files for OpenCV 2.4.8 in the project. The other OpenCV versions work fine too, but you need to change the OpenCV setting in the project.
   b. Modify the configurations for sequences in "example.cpp" £¨we give an example for the sequence "torus"£©.
   c. Set the initial bounding box for target in "init.txt" and generate the image list in "imglist.txt". The initialization bounding box is in the formate of (x-topleft,y-topleft,width,height), followed by an "enter" for data reading. 
   d. Complie and run.
   
   Note: the file "PARAMETERS.txt" contains some parameter settings we use in our experiments. Since we have changed some parts of the codes, you may not get the exact same results we have in [1]. But these settings work well too.
   
   This software is tested under 64bit windows and Visual Studio 2012 compiler. It should be used only for research purposes. 
   For any questions/remarks/problems with the code please contact zwcai@ucsd.edu or lywen.cv.workbox@gmail.com.

==============================================
UPDATE LOG
==============================================
v2.1 2015.5.15
Some useless files and variables are removed for better readability. 
Memory release problem is solved to some degree. But the memory usage is still not controlled very well. The major problem is within online SVM files. 
Parameter setting is introduced for setting parameters on other untested sequences.  

v2.0 2015.4.18
More comments are added.
The names of some variable are modified to improve the readability.
The online SVM model is re-initialized every several frames because of runing speed.

v1.0 2015.4.16
Preliminary version released.

##################################################################


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


##################################################################