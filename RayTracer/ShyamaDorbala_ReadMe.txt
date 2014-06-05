CSCI 420
Assignment #3: Ray tracing

FULL NAME: Shyama Dorbala
USC ID: 6120-0631-99


USAGE
------------------
1. For seeing still ray traced images
	 make clean 
	 make
	 ./ShyamaDorbala_assign3 <scenefile>

2. For seeing animation based on scene changes
	 make clean
	 make
	 ./ShyamaDorbala_assign3 <scenefile> 0 1 0 0
	 **P.S. For seeing this properly, its recommended to use spheres.scene

3. For seeing Motion Blur
	 make clean
	 make
	 ./ShyamaDorbala_assign3 <scenefile> 1 0 0 0
	 **P.S. For seeing this motion, its recommended to use table.scene

4. For soft shadows
	 make clean
	 make
	 ./ShyamaDorbala_assign3 <scenefile> 0 0 1 0
	 **P.S. For seeing this motion, its recommended to use test2.scene


MANDATORY FEATURES
------------------

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  Yes

2) Ray tracing sphere                     Yes

3) Triangle Phong Shading                 Yes

4) Sphere Phong Shading                   Yes

5) Shadows rays                           Yes

6) Still images                           Yes
   
7) Extra Credit

		a. Animations by changing the scene
			 - The lights in the scene are constantly moved 
				 to create animations
			 
		b. Motion Blur animations
			 - This is implemented using ray tracing and gl accumulation buffers
			 - The accumulation buffer accumulates the color for a set of frames
				 and then is returned to the display.
			 - This gives a motion blur effect
		
		c. Soft shadows
			 - Used an area light source, to shoot different shadow rays and
				 average the values to form soft shadows.


SUBMITTED IMAGES
------------------

Different kinds of images are submitted:
	1. Still images with ray tracing for various scenes
		 	- with fov as 60 and 90
			- with soft shadows, as the number of components to make 
				area light increases
	2. light motion animation frames: lM_001.jpg..
	3. motion blur animation frames: mB_001.jpg..
			
TIMING
------------------

For the still images, here are the details of the computation times:
	1. test2.scene: 3 sec
	2. table.scene: 45 sec - 1 min
	3. siggraph.scene: 6 min - 7 min
