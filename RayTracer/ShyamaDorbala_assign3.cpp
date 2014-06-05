/*
CSCI 480
Assignment 3 Raytracer
Name: Shyama Dorbala
USC ID: 6120-0631-99
*/

/* Remove the warnings from parse, etc. */
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wunused-result"

/* Standard CPP Libraries */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <string.h>
#include <sstream>
#include <vector>
using namespace std;

/* GL Libraries */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <pic.h>

/* Global Defines */
#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 100

/* Filename handle */
char *filename=0;

/* Different display modes */
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

/* Size of the window */
#define WIDTH 640
#define HEIGHT 480

/* Field of view of the camera */
#define fov 90.0
#define PI 3.14159265

/* Enable Motion Blur */
bool motionBlur=0;

/* Enable light motion for animation */
bool lightMotion=0;

/* Enable soft shadows */
bool softShadow=0;
int num_areaLightComps=35;

/* Max co-ords in world space based on fov */
double yMax=tan((double) PI*fov/(2*180));
double xMax=yMax*((double)WIDTH)/((double)HEIGHT);

/* Color buffer to hold the rgb values each 8-bit */
unsigned char buffer[HEIGHT][WIDTH][3];

/* Program-specific structs for easier access */
struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct point
{
	double x;
	double y;
	double z;
};

/* In-out point struct */
struct inOut
{
	bool in; // =0 if outside the triangle, =1 if inside
	double bary[3]; // bary centric co-ordinates for a given point
};

/* Intersection point struct */
struct intxPoint
{
	point p; // intersection point
	double t; // t value of the intersection point
	int tID; // index of the object in the array of objects (object=triangle/sphere)
	int tObj; // =1 if the point is on sphere, = 2 if triangle
	inOut iO; // contains the barycentric co-ordinates for the intxPoint
};

/* Objects */
typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

/* Global variables */
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];

double ambient_light[3];

struct point cam;

int num_triangles=0;
int num_spheres=0;
int num_lights=0;


/* pixel handling function prototypes */
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

/* Find difference of two point types */
point findDiff(point A,point B)
{
	point C;

	C.x=A.x-B.x;
	C.y=A.y-B.y;
	C.z=A.z-B.z;

	return C;
}

/* Divide a point type by a scalar */
point scalarDiv(point A,double a)
{
	point B;

	B.x=0.0;
	B.y=0.0;
	B.z=0.0;

	if (abs(a)>1e-10) //a>0
	{
		B.x=A.x/a;
		B.y=A.y/a;
		B.z=A.z/a;
	}

	return B;
}

/* Find magnitude of a point */
double findMag(point A)
{
	double Amag;

	Amag = sqrt(pow(A.x,2)+pow(A.y,2)+pow(A.z,2));
	return Amag;
}

/* Normalize a point */
point unitize(point A)
{
	point U;
	double Amag;
	
	Amag=findMag(A);
	U = scalarDiv(A,Amag);

	return U;
}

/* Find cross-product of two vectors A,B */
point crossProduct(point A,point B)
{
	point C;

	C.x = (A.y*B.z - B.y*A.z);
	C.y = (B.x*A.z - A.x*B.z);
	C.z = (A.x*B.y - A.y*B.x);
	
	return C;
}

/* Find dot-product of two vectors A,B */
double dotProduct(point A,point B)
{
	double C;

	C = (A.x*B.x + A.y*B.y + A.z*B.z);
	return C;
}

/* Find a point on the ray that originates from src, 
	 directed along dir with parameter t 
*/
point shootRay(point src,point dir,double t)
{
	point p;
	
	p.x=src.x+t*(dir.x);
	p.y=src.y+t*(dir.y);
	p.z=src.z+t*(dir.z);

	return p;
}

/* Find the parameter t on the ray originating from src,
	 directed along dir, which corresponds to the point 
	 that intersects the given sphere
*/
double intersectSphere(Sphere sphere,point src,point dir)
{
	double b,c,t,t1,t2;
	t1=0;
	t2=0;

	// From the ray intersects sphere equation	
	c = pow((src.x-sphere.position[0]),2) 
			+ pow((src.y-sphere.position[1]),2) 
			+ pow((src.z-sphere.position[2]),2) 
			- pow(sphere.radius,2);
	b = 2*(dir.x*(src.x-sphere.position[0]) 
			+ dir.y*(src.y-sphere.position[1]) 
			+ dir.z*(src.z-sphere.position[2]));

	// Check if the determinant is positive	
	if ( (pow(b,2)-4*c)>0 ) 
	{
		
		// Points of intersection
		t1=(((-1)*b)+sqrt(pow(b,2)-4*c))/2;
		t2=(((-1)*b)-sqrt(pow(b,2)-4*c))/2;
		
		// Min positive t
		t = min( t1, t2 );
		
		// Ignore if negative	
		if (t<0) 
			t=-1;

		// If one of the positive t values is close to 0, take the other value as a valid t
		else if (t<1e-10)
		{
			if (t1<1e-15 && t2>1e-15) 
				t=t2;
			else if (t2<1e-15 && t1>1e-15) 
				t=t1;
			else t=-1; // ignore if both are close to 0
		}
	}
	else t=-1; // ignore if determinant is negative

	return t;
}

/* Find the line segment formed by vertices V1 and V2 */
point getSide(Vertex v1, Vertex v2)
{
	point AB;
	
	AB.x=v1.position[0]-v2.position[0];
	AB.y=v1.position[1]-v2.position[1];
	AB.z=v1.position[2]-v2.position[2];
	
	return AB;
}

/* Check if two vectors A and B are equal or same */
bool checkEqual(point A,point B)
{
	bool eq;
	
	if ((abs(A.x-B.x)<1e-10) && (abs(A.y-B.y)<1e-10) && (abs(A.z-B.z)<1e-10))
		eq=1;
	else eq=0;

	return eq;
}

/* Check if a point p lies inside the given triangle */
inOut inoutTest(Triangle triangle,point p)
{
	inOut iO;
	Vertex P;
	point PA,PB,PC;
	point uPAxPB,uPBxPC,uPCxPA,PAxPB,PBxPC,PCxPA;
	double A,B,C;	
	
	// Convert point to a vertex	
	P.position[0]=p.x;
	P.position[1]=p.y;
	P.position[2]=p.z;
	
	// Find the sides of the triangle
	PA=getSide(P,triangle.v[0]);
	PB=getSide(P,triangle.v[1]);
	PC=getSide(P,triangle.v[2]);

	// Find the un-normalized cross products of the sides
	uPAxPB=crossProduct(PA,PB);
	uPBxPC=crossProduct(PB,PC);
	uPCxPA=crossProduct(PC,PA);
	
	// Magnitudes of the cross products
	C=findMag(uPAxPB);
	A=findMag(uPBxPC);
	B=findMag(uPCxPA);	

	// Find the barycentric co-ords
	if ((A+B+C)>1e-20)	
	{
		iO.bary[0]=A/(A+B+C);//alpha
		iO.bary[1]=B/(A+B+C);//beta
		iO.bary[2]=C/(A+B+C);//gamma
	} 

	// Normalize the cross products	
	PAxPB=scalarDiv(uPAxPB,C);	
	PBxPC=scalarDiv(uPBxPC,A);	
	PCxPA=scalarDiv(uPCxPA,B);	

	// Find if the normalized cross products are equal, since they are unit vectors now
	if(checkEqual(PAxPB,PBxPC) && checkEqual(PAxPB,PCxPA))
		iO.in=1;
	else iO.in=0;

	return iO;				
}

/* Find the parameter t corresponding to a point on the ray which originates 
	 from src, along dir and intersects the triangle 
*/
double intersectTriangle(Triangle triangle,point src,point dir)
{
	point AB,AC,n;	
	double d,t;
	t=0;
	d=0;
	
	// Get sides of the triangle to find the normal	
	AB=getSide(triangle.v[0],triangle.v[1]);
	AC=getSide(triangle.v[0],triangle.v[2]);

	// Cross Product to find normal	
	n=crossProduct(AB,AC);
	n=unitize(n);

	// d of the ray intersects plane equation
	d = (n.x)*(-1)*triangle.v[1].position[0]
			- (n.y)*triangle.v[1].position[1]
			- (n.z)*triangle.v[1].position[2];

	// If denominator is close to 0, ignore
	if (abs(dotProduct(n,dir))<1e-35) 
		t=-1;
	else
		t=(-1)*(dotProduct(n,src)+d)/(dotProduct(n,dir));
	
	return t;
}

/* Find the normal of the sphere with index tID at the point p */
point findSphereNormal(point p,int tID)
{
	point n;

	// based on the equation
	n.x=(p.x-spheres[tID].position[0])/spheres[tID].radius;
	n.y=(p.y-spheres[tID].position[1])/spheres[tID].radius;
	n.z=(p.z-spheres[tID].position[2])/spheres[tID].radius;
	return n;
}

/* Interpolate various entities of a triangle based on 
	 the barycentric co-ordinates given by the inOut object,
	 ID = 0, interpolate the normal
			= 1, interpolate the diffuse component,
			= 2, interpolate the specular component 
*/
point interpolate(Triangle triangle,inOut iO,int ID)
{
	point P;
	
	if (ID==0)
	{
		P.x = iO.bary[0]*triangle.v[0].normal[0]
				 +iO.bary[1]*triangle.v[1].normal[0]
				 +iO.bary[2]*triangle.v[2].normal[0];

		P.y = iO.bary[0]*triangle.v[0].normal[1]
				 +iO.bary[1]*triangle.v[1].normal[1]
				 +iO.bary[2]*triangle.v[2].normal[1];

		P.z = iO.bary[0]*triangle.v[0].normal[2]
				 +iO.bary[1]*triangle.v[1].normal[2]
				 +iO.bary[2]*triangle.v[2].normal[2];
	}

	else if (ID==1)
	{
		P.x = iO.bary[0]*triangle.v[0].color_diffuse[0]
				 +iO.bary[1]*triangle.v[1].color_diffuse[0]
				 +iO.bary[2]*triangle.v[2].color_diffuse[0];

		P.y = iO.bary[0]*triangle.v[0].color_diffuse[1]
				 +iO.bary[1]*triangle.v[1].color_diffuse[1]
				 +iO.bary[2]*triangle.v[2].color_diffuse[1];

		P.z = iO.bary[0]*triangle.v[0].color_diffuse[2]
				 +iO.bary[1]*triangle.v[1].color_diffuse[2]
				 +iO.bary[2]*triangle.v[2].color_diffuse[2];
	}

	else if (ID==2)
	{
		P.x = iO.bary[0]*triangle.v[0].color_specular[0]
				 +iO.bary[1]*triangle.v[1].color_specular[0]
				 +iO.bary[2]*triangle.v[2].color_specular[0];

		P.y = iO.bary[0]*triangle.v[0].color_specular[1]
				 +iO.bary[1]*triangle.v[1].color_specular[1]
				 +iO.bary[2]*triangle.v[2].color_specular[1];

		P.z = iO.bary[0]*triangle.v[0].color_specular[2]
				 +iO.bary[1]*triangle.v[1].color_specular[2]
				 +iO.bary[2]*triangle.v[2].color_specular[2];
	}

	return P;
}

/* Compute phong shading color at point p, on an object identified by tObj,
	 indexed by tID, lit by light and viewed from cam
*/
point phong(point p,int tID,int tObj,inOut iO,Light light,point camera)
{
	point n,l,v,r,kd,ks;
	point I; 
	double lDotN,rDotV,alpha;

	// Convert to a point type
	l.x=light.position[0];
	l.y=light.position[1];
	l.z=light.position[2];

	// Find the direction vector from given point to the given light
	l=unitize(findDiff(l,p));
	
	// Find the direction vector from the given point to the given cam position
	v=unitize(findDiff(camera,p));

	// Find the normal, kd, ks, shininess based on which object the point intersects
	if (tObj==1) // if sphere
	{
		n=findSphereNormal(p,tID);
		
		kd.x=spheres[tID].color_diffuse[0];	
		kd.y=spheres[tID].color_diffuse[1];	
		kd.z=spheres[tID].color_diffuse[2];	
		
		ks.x=spheres[tID].color_specular[0];	
		ks.y=spheres[tID].color_specular[1];	
		ks.z=spheres[tID].color_specular[2];
		
		alpha=spheres[tID].shininess;	
	}

	else if (tObj==2) // if triangle
	{
		n=unitize(interpolate(triangles[tID],iO,0));
		kd=interpolate(triangles[tID],iO,1);
		ks=interpolate(triangles[tID],iO,2);
		
		alpha=iO.bary[0]*triangles[tID].v[0].shininess
				 +iO.bary[1]*triangles[tID].v[1].shininess
				 +iO.bary[2]*triangles[tID].v[2].shininess;
	}

	// Compute dot product l.n, clamp it to 0-1	
	lDotN=dotProduct(l,n);
	if (lDotN<0)
		lDotN=0;	
	else if (lDotN>1.f) 
		lDotN=1.f;
	
	// Find the reflection vector with l and n as r=2(l.n)n-l
	r.x=2*lDotN*n.x-l.x;
	r.y=2*lDotN*n.y-l.y;
	r.z=2*lDotN*n.z-l.z;

	// Compute dot product r.v, clamp it to 0-1	
	rDotV=dotProduct(r,v);
	if (rDotV<0) 
		rDotV=0;	
	else if (rDotV>1.f) 
		rDotV=1.f;

	// Compute Intensity using the Phong equation
	I.x=light.color[0]*((kd.x)*lDotN+((ks.x)*pow((rDotV),(alpha)))); // r
	I.y=light.color[1]*((kd.y)*lDotN+((ks.y)*pow((rDotV),(alpha)))); // g 
	I.z=light.color[2]*((kd.z)*lDotN+((ks.z)*pow((rDotV),(alpha)))); // b
	return I;
}

/* Find the point of intersection of the rays originating from pNear, 
	 going to pFar.
	 If shadow = 0, pFar=pixel location, pNear=camera, ray's origin = pFar
						 = 1, pFar=light, pNear=point of intersection, ray's origin = pNear
	 intxPoint gives the point of intersection, the object of intersection, 
	 and the parameter t corresponding to the point on the ray, 
	 and if the object is a triangle, then the corresponding barycentric co-ords
*/
intxPoint intersectObjects(point pFar,point pNear,int shadow)
{
	point p,q,raySrc,dir,pixPoint;
	intxPoint intxObj;
	inOut iO;	
	double t,tMin,tMax,tS,tT;
	int tID,tObj;
	
	// Initialize t, index and object
	tMin=0;
	tID=-1;
	tObj=-1;	

	// Change origin of the ray based on if its a shadow ray or not	
	if (shadow==0)		
		raySrc=pFar; // pNear is cam
	else if (shadow==1)
		raySrc=pNear; // pFar is light
	
	// Find the direction vector of the ray	
	dir.x=pFar.x-pNear.x;
	dir.y=pFar.y-pNear.y;
	dir.z=pFar.z-pNear.z;
	dir=unitize(dir);

	// Find the point of intersection with objects of the ray starting from raySrc along dir
	// If it is a shadow ray, the point neednt be found, just the t value would suffice 
	// because there is no phong shading for points in shadow
	
	// First, with spheres	
	for (int i=0;i<num_spheres;i++)
	{
		tS=intersectSphere(spheres[i],raySrc,dir); // Find t
		if (tMin==0 && tS>1e-10) // Check for positive t 
		{
			tMin=tS;
			tID=i;
			tObj=1;
		}
		else if (tS<=tMin && tS>1e-10) // Minimum positive t of all spheres 
		{
			tMin=tS;
			tID=i;
			tObj=1;
		}
	}
	// with Triangles
	for (int i=0;i<num_triangles;i++)
	{
		tT=intersectTriangle(triangles[i],raySrc,dir); // Find t
		p=shootRay(raySrc,dir,tT); // Find the point corresponding to tT
		iO=inoutTest(triangles[i],p); // Check if it lies in the triangle
		
		if (iO.in==1) // If inside the triangle
		{	
			if (tMin==0 && tT>1e-5) // Check for positive t
			{
				tMin=tT;
				tID=i;
				tObj=2;
				if (shadow==0)
					q=p;
			  intxObj.iO.bary[0]=iO.bary[0];	
			  intxObj.iO.bary[1]=iO.bary[1];	
			  intxObj.iO.bary[2]=iO.bary[2];	
			}
			else if (tT<tMin && tT>1e-5) // Minimum positive t of all triangles and spheres 
			{
				tMin=tT;
				tID=i;
				tObj=2;
				if (shadow==0)
					q=p;
			  intxObj.iO.bary[0]=iO.bary[0];	
			  intxObj.iO.bary[1]=iO.bary[1];	
			  intxObj.iO.bary[2]=iO.bary[2];	
			}
		}
	}	

	// Check if the intersection point lies beyond the light in case of a shadow ray
	if (shadow==1) 
	{
		// Find t on the ray that corresponds to the light position
		if (dir.x!=0)
		{
			tMax=(pFar.x-raySrc.x)/dir.x;
		}
		else if (dir.y!=0)
		{
			tMax=(pFar.y-raySrc.y)/dir.y;
		}
		else if (dir.x!=0)
		{
			tMax=(pFar.z-raySrc.z)/dir.z;
		}
		else tMax=0;
	
		// The point should lie on the ray before the light
		// thus, t of the Point should be less than the t of the light	
		if (tMin>=tMax)
		{
			tObj=-1;
			tID=-1;
		}
	}
	// If its not a shadow ray and the intersection is with circle
	// find the point of intersection, wasnt found so far for circle
	// If its a triangle, the point is already found
	else if ((tMin>=0) && (tObj==1))
		q=shootRay(raySrc,dir,tMin); 
	
	intxObj.p=q;
	intxObj.t=tMin;
	intxObj.tID=tID;
	intxObj.tObj=tObj;
	return intxObj;
}

/* Find the color of the pixel at (x,y) in the image space */
point findColor(int x, int y) 
{
	point p,q,dir,light,lightS;
	point black,pixColor,temp,tempN;
	intxPoint intxObj,intxShadow, intxSoft;
	int hCtr=1;
	
	// Counter to loop through lights
	if (softShadow)
		hCtr = num_areaLightComps+1;	

	// Default color	
	black.x=0.0;
	black.y=0.0;
	black.z=0.0;
	pixColor=black;

	// Convert pixel co-ordinates to world co-ordinates	
	p.x=(((double)x/(double)WIDTH)*2*xMax)-xMax;
	p.y=(((double)y/(double)HEIGHT)*2*yMax)-yMax;
	p.z=-1;	

	// Shoot rays from pixel point p to intersect objects
	intxObj=intersectObjects(p,cam,0);
	
	// If it intersects any object only, enter the loop to find phong color
	if (intxObj.tID!=-1) 
	{	
		// Point of intersection
		q=intxObj.p; 
		
		// Add the global ambient light first
		pixColor.x+=ambient_light[0];
		pixColor.y+=ambient_light[1];
		pixColor.z+=ambient_light[2];
		
		// Find the contribution of each light source
		for (int h=0;h<num_lights;h+=hCtr)
		{
			// Convert to a point type
			light.x=lights[h].position[0];
			light.y=lights[h].position[1];  
			light.z=lights[h].position[2];
		
			// Shoot shadow rays to the light from point and find intersection with the objects
			intxShadow=intersectObjects(light,q,1);
		
			// If it doesnt intersect any object only, enter to find contribution of the light
			if ((intxShadow.tID==-1))
			{
				// Soft Shadow Computation: shoot rays around the light and average
				if (softShadow) 
				{
					temp=black;
					
					for (int j=h;j<(h+(num_areaLightComps)+1);j++)
					{
						// Convert to a point type
						lightS.x=lights[j].position[0];
						lightS.y=lights[j].position[1];  
						lightS.z=lights[j].position[2];
					
						// Find intersection with Area light	
						intxSoft=intersectObjects(lightS,q,1);
					
						if (intxSoft.tID==-1)
						{
							// Calculate phong color for each of the lights
							tempN=phong(q,intxObj.tID,intxObj.tObj,intxObj.iO,lights[j],cam);
							
							// Sum the results
							temp.x+=tempN.x;
							temp.y+=tempN.y;
							temp.z+=tempN.z;
						}
					}
					
					// Average the area light
					if (num_areaLightComps!=0)
						temp=scalarDiv(temp,num_areaLightComps+1);
				}				
				else
					// Find the phong model color
					temp=phong(q,intxObj.tID,intxObj.tObj,intxObj.iO,lights[h],cam);
				
				// Add to the pixel color
				pixColor.x+=temp.x;
				pixColor.y+=temp.y;
				pixColor.z+=temp.z;
			}
		}

		// If pixel color goes above 1.f, clamp it to 1.f
		if (pixColor.x > 1) pixColor.x=1.f;
		if (pixColor.y > 1) pixColor.y=1.f;
		if (pixColor.z > 1) pixColor.z=1.f;
	}
	else pixColor=black;
	
	return pixColor;
}

/* Fill the color buffer with all the computed values for a particular scene */
void fill_color()
{
  unsigned int x,y;
	point pixColor;
 
	// For each pixel, find the color 
	for(x=0; x<WIDTH; x++)
  {
    for(y=0;y < HEIGHT;y++)
    {
			// Find the color
      pixColor=findColor(x,y);

			// Put the color values in the buffer
			plot_pixel_jpeg(x,y,abs(pixColor.x)*255,abs(pixColor.y)*255,abs(pixColor.z)*255);
    }
  }
}

/* With the values in the color buffer, draw them to the GL window */
void draw_scene()
{
  unsigned int x,y;
  glPointSize(2.0);
  
  glBegin(GL_POINTS);  
	for(x=0; x<WIDTH; x++)
  {
    for(y=0;y < HEIGHT;y++)
    {
			plot_pixel_display(x,y,buffer[HEIGHT-y-1][x][0],buffer[HEIGHT-y-1][x][1],buffer[HEIGHT-y-1][x][2]);
    }
  }
  glEnd();

	if(!motionBlur && !lightMotion)
 		glFlush(); // Comment this
  printf("Done!\n"); 
	if(!motionBlur && !lightMotion)
		fflush(stdout); // Comment this
}

/* Plot to the GL window */
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

/* Put the values in a buffer */
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

/* Choose to draw to window or save to jpeg */
void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

/* save the color buffer to jpeg */
void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

	if (motionBlur)
	{
		for (int i=HEIGHT-1; i>=0; i--) 
		{
    	glReadPixels(0, HEIGHT-1-i, WIDTH, 1, GL_RGB, GL_UNSIGNED_BYTE,
                 &in->pix[i*in->nx*in->bpp]);
  	}
	}
	else
		memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  
	if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

/* Reading the scene files and error checking */
void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l,l1;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
		
		/* adding some extra lights to make it an area light */
		l1=l;
		if (softShadow)
		{
			for (int aL=1;aL<(num_areaLightComps+1);aL++)
			{
				l1.position[2]+=(aL*0.001-0.02);
	  		lights[num_lights++] = l1;
			}
		}
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void idle()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_ACCUM_BUFFER_BIT);
}

/* display to the window 
	 Three kinds of displays:
	 1. when motionBlur=0 and lightMotion=0 : A still image is displayed after ray tracing
	 2. when motionBlur=0 and lightMotion=1 : An animation is produced as the lights in the scene move
	 3. when motionBlur=1 and lightMotion=0 : An animation is produced showcasing the motion blur effect
*/
void display()
{
  static int once=0;
  static int count=0;
  static int n=15;
  static double blur=0;
	static int k=0;
	string fnameMB;	
	string fname;	

	// If light is moving, every frame is refreshed
	if (lightMotion)
		glClear(GL_COLOR_BUFFER_BIT);// | GL_DEPTH_BUFFER_BIT);

	// Camera Position for the scene
	cam.x=0.0;
 	cam.y=0.0;
 	cam.z=0.0;

	// Reset the world
	glLoadIdentity();

	// First time writing to color Buffer and lightMotion computation
	if((lightMotion && once<100) || (!once && !lightMotion))
  {
		// First time filling the color buffer by ray Tracing	
		fill_color();
		
		// If motionBlur=1, dont draw it here
		if (!motionBlur)
			draw_scene();

		// Save lightMotion screenshots
		if (mode==MODE_JPEG & !motionBlur & lightMotion)
		{
			k++;
			
			if (k>0 && k<100) 
			{
				stringstream ss;
    			ss << k;
    			ss >> fnameMB;
	
				if (k<10)
					fnameMB.insert(0,"00");
				else if (k<100)
					fnameMB.insert(0,"0");
				
				fname="lM_"+fnameMB+".jpg";
				
				strcpy(filename, fname.c_str());
				
				save_jpg();
			}
		}
		else if (mode==MODE_JPEG && !motionBlur && !lightMotion)
			save_jpg();
		
		// Move the lights in the scene by an offset
		if (lightMotion && !motionBlur)
		{
			for (int l=0;l<num_lights;l++)
			{
				lights[l].position[0]+=(once*0.1-2);
				lights[l].position[0]+=(once*0.1-2);
			}
		}
  }
  once++;
	
	// Swap buffers and Redisplay for lightMotion	
	if (lightMotion && once<100) 
	{
		glutSwapBuffers();
		glutPostRedisplay();
	}		

	// Motion Blur computation
	if (motionBlur)
	{
		// variable to keep a track of iterations
		count++;

		// Measure of the movement/blur
		blur+=2;
		
		// Move the scene
		glTranslatef(blur,0,0);
		
		// Redraw the scene
		draw_scene();

		// Load into the accumulation buffer and keep accumulating a bunch of frames
		if (once==1)
			glAccum(GL_LOAD,0.5/n); //0.5/n
		else 
			glAccum(GL_ACCUM,0.5/n);

		// Stop the blurring when count reaches 2*n
		if (count<2*n) 
		{ 
			// Copy the accumulation buffer contents
			glAccum(GL_RETURN,1.0f);

			// Swap buffers and Redisplay
			glutSwapBuffers();
			glutPostRedisplay();

			// Save motionBlur images
			if (mode==MODE_JPEG)
			{
				k++;

				if (k>0 && k<100) 
				{
					stringstream ss;
    				ss << k;
    				ss >> fnameMB;
	
					if (k<10)
						fnameMB.insert(0,"00");
					else if (k<100)
						fnameMB.insert(0,"0");
					
					fname="mB_"+fnameMB;
					fname+= ".jpg";
	
					strcpy(filename, fname.c_str());

					save_jpg();
				}
				glClear(GL_COLOR_BUFFER_BIT);
			}
		}

		// Bunch of frames to accumulate, the accumulation buffer is reloaded
		// when once=1
		if (once>=n) {
			once=1;
		}
	}

}

int main (int argc, char ** argv)
{
	
  if (argc!=2 && argc!=3 && argc!=6)
  {  
    // For animation, MUST give all flags
		printf ("\nUsage1: %s <scenefile> <motionBlurFlag> <lightMotionFlag> <softShadowFlag> <jpegFlag>\n", argv[0]);
		printf ("\nFlag=1, enable\nFlag=0, disable");
		printf ("\nPlease enable only one of motion Blur, light motion or soft shadows at one time\n");
		printf ("\neg for motion blur: ./ShyamaDorbala_assign3 table.scene 0 1 0 0");
		printf ("\neg for light motion: ./ShyamaDorbala_assign3 spheres.scene 1 0 0 0");
		printf ("\neg for soft shadows: ./ShyamaDorbala_assign3 test2.scene 0 0 1 0\n\n");
		// For simple ray tracing
    printf ("\nUsage2: %s <scenefile> [jpegname]\n\n", argv[0]);
    exit(0);
  }

  if(argc == 6)
  {
  	if ( atoi(argv[4])==1 )
			mode = MODE_JPEG;
		else 
			mode = MODE_DISPLAY;
		motionBlur = atoi(argv[2]);
		lightMotion = atoi(argv[3]);
		softShadow = atoi(argv[4]);
		filename = argv[5];
  }
	else if (argc == 3)
	{
    mode = MODE_JPEG;
    filename = argv[2];
	}
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

	// Double buffering enabled in case of animations	
	if (motionBlur || lightMotion)
  	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ACCUM); 
	else
  	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE | GLUT_ACCUM);
 
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Shyama's Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
