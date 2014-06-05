/*
  CSCI 480
  Assignment 2
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <string.h>
#include <sstream>
#include <vector>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "pic.h"

using namespace std;

GLint g_vMousePos[2] = {0, 0};
GLint vMouseDelta[2] = {0, 0};
GLint g_iLeftMouseButton = 0;    
GLint g_iMiddleMouseButton = 0;
GLint g_iRightMouseButton = 0;
GLint g_iMenuId;

int a=0;
int b=0;
vector<double> u;
static char * fn = new char[8];

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;
CONTROLSTATE g_ControlState = ROTATE;

/* state of the world */
GLfloat g_vLandRotate[3] = {0.0, 0.0, 0.0};
GLfloat g_vLandTranslate[3] = {0.0, 0.0, 0.0};
GLfloat g_vLandScale[3] = {1.0, 1.0, 1.0};
GLfloat g_incTranslate=0;
GLfloat g_keyTranslate[3]={0.0,0.0,0.0};
GLint g_pause=0;
GLint g_MoveCam=1;
GLint g_Move=0;
GLint g_Sphere=0;
GLint g_NoSphere=0;
GLint g_takeShot=0;

Pic * texData;

/* represents one control point along the spline */
struct point {
   double x;
   double y;
   double z;
};

/* spline struct which contains how many control points, and an array of control points */
struct spline {
   int numControlPoints;
   struct point *points;
};

/* the spline array */
struct spline *g_Splines;
struct spline *g_subSplines;
struct spline *g_subTangents;
struct spline *g_subNormals;
struct spline *g_subBiNormals;

/* total number of splines */
int g_iNumOfSplines;
int g_iNumOfSubSplines;
double dist = 0.01;
double maxHeight;
double maxLineLength=0.001;

struct skyCube {
	GLuint skyTexture[6];
	GLubyte * sky_texels[6];
	GLuint tWidth[6];
	GLuint tHeight[6];
};

struct skyCube my_sky;

struct crossBar {
	GLuint crossTexture;
	GLubyte * cross_texels;
	GLuint crossWidth;
	GLuint crossHeight;
};

struct crossBar my_bar;
struct crossBar my_toy;

double max(double a,double b,double c)
{
	if(!c)
	{
		if(a>b) return a;
		else return b;
	}
	else if (a>b) return max(a,c);
	else return max(b,c);
}

double find_maxCoord(spline g_aSpline)
{
	double maxCord=0;
	double maxPoint=0;
	point pt;	
	
	for(int i=0;i<g_aSpline.numControlPoints;i++)
	{
		if(g_aSpline.points[i].x<0) 
			pt.x=-g_aSpline.points[i].x; 
		else 
			pt.x=g_aSpline.points[i].x;
		
		if(g_aSpline.points[i].y<0) 
			pt.y=-g_aSpline.points[i].y; 
		else 
			pt.y=g_aSpline.points[i].y;
		
		if(g_aSpline.points[i].z<0) 
			pt.z=-g_aSpline.points[i].z; 
		else 
			pt.z=g_aSpline.points[i].z;
		
		maxPoint= max(pt.x,pt.y,pt.z);
		maxCord= max(maxPoint,maxCord);
	}
	return maxCord;
}

/* Write a screenshot to the specified filename */
void saveScreenshot (char *filename)
{
  int i, j;
  Pic *in = NULL;

  if (filename == NULL)
    return;

  /* Allocate a picture buffer */
  in = pic_alloc(640, 480, 3, NULL);

  printf("File to save to: %s\n", filename);

  for (i=479; i>=0; i--) {
    glReadPixels(0, 479-i, 640, 1, GL_RGB, GL_UNSIGNED_BYTE,
                 &in->pix[i*in->nx*in->bpp]);
  }

  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}

/* loads splines from track.txt */
int loadSplines(char *argv) {
  char *cName = (char *)malloc(128 * sizeof(char));
  FILE *fileList;
  FILE *fileSpline;
  int iType, i = 0, j, iLength;
	double maxCord;

  /* load the track file */
  fileList = fopen(argv, "r");
  if (fileList == NULL) {
    printf ("can't open file\n");
		exit(1);
  }
  
  /* stores the number of splines in a global variable */
  int a1=fscanf(fileList, "%d", &g_iNumOfSplines); //FIXME: variable name alert!

  g_Splines = (struct spline *)malloc(g_iNumOfSplines * sizeof(struct spline)); 
	
	/* reads through the spline files */
  for (j = 0; j < g_iNumOfSplines; j++) {
    i = 0;
    a1=fscanf(fileList, "%s", cName);
    fileSpline = fopen(cName, "r");
		cout << "my file " << cName << "\n";

    if (fileSpline == NULL) {
      printf ("can't open file in cName\n");
			cout << "my file in NULL" << cName << "\n";
			exit(1);
    }

    /* gets length for spline file */
    a1=fscanf(fileSpline, "%d %d", &iLength, &iType);

    /* allocate memory for all the points */
    g_Splines[j].points = (struct point *)malloc(iLength * sizeof(struct point));
    g_Splines[j].numControlPoints = iLength;


		/* saves the data to the struct */
    while (fscanf(fileSpline, "%lf %lf %lf", 
	   &g_Splines[j].points[i].x, 
	   &g_Splines[j].points[i].y, 
	   &g_Splines[j].points[i].z) != EOF) {
      i++;
    }
		
		fclose(fileSpline);
		maxCord=find_maxCoord(g_Splines[j]);
		for(int q=0;q<g_Splines[j].numControlPoints;q++)
		{
			g_Splines[j].points[q].x = g_Splines[j].points[q].x/maxCord;
			g_Splines[j].points[q].y = g_Splines[j].points[q].y/maxCord;
			g_Splines[j].points[q].z = g_Splines[j].points[q].z/maxCord;
		}
  }
	
  free(cName);

  return 0;
}

double find_maxHeight()
{
	double maxHeight=0;
	for (int j = 0; j < g_iNumOfSubSplines; j++) 
	{
		for(int q=0;q<g_subSplines[j].numControlPoints;q++)
		{
			maxHeight=max(maxHeight,g_subSplines[j].points[q].y);
		}
  }
	return maxHeight;
}

/* finds the interpolated point for a particular u 
 * between points p1 and p2, with tangents p3 and p4
 */
point find_point(double u, point p1, point p2, point p3, point p4)
{
	point p;
	
	p.x =	(double) pow(u,3)*(-0.5*p1.x + 1.5*p2.x - 1.5*p3.x + 0.5*p4.x) + 
				pow(u,2)*(p1.x - 2.5*p2.x + 2.0*p3.x - 0.5*p4.x) + 
				u*(-0.5*p1.x + 0.5*p3.x) + (p2.x);
	
	p.y = (double) pow(u,3)*(-0.5*p1.y + 1.5*p2.y - 1.5*p3.y + 0.5*p4.y) + 
				pow(u,2)*(p1.y - 2.5*p2.y + 2.0*p3.y - 0.5*p4.y) + 
				u*(-0.5*p1.y + 0.5*p3.y) + (p2.y);
	
	p.z = (double) pow(u,3)*(-0.5*p1.z + 1.5*p2.z - 1.5*p3.z + 0.5*p4.z) + 
				pow(u,2)*(p1.z - 2.5*p2.z + 2.0*p3.z - 0.5*p4.z) + 
				u*(-0.5*p1.z + 0.5*p3.z) + (p2.z);
	
	return p;
}

point unitize(point A)
{
	point U;
	double Amag;
	
	Amag = sqrt(pow(A.x,2)+pow(A.y,2)+pow(A.z,2));
	if (Amag>0)
	{
		U.x=A.x/Amag;
		U.y=A.y/Amag;
		U.z=A.z/Amag;
	}
	else
	{
		U.x=0;
		U.y=0;
		U.z=0;
	}

	return U;
}

point find_tangent(double u, point p1, point p2, point p3, point p4)
{
	point p;
	
	p.x = 3*pow(u,2)*(-0.5*p1.x + 1.5*p2.x - 1.5*p3.x + 0.5*p4.x) + 
				2*u*(p1.x - 2.5*p2.x + 2.0*p3.x - 0.5*p4.x) + 
				(-0.5*p1.x + 0.5*p3.x);// + (p2.x);
	
	p.y = 3*pow(u,2)*(-0.5*p1.y + 1.5*p2.y - 1.5*p3.y + 0.5*p4.y) + 
				2*u*(p1.y - 2.5*p2.y + 2.0*p3.y - 0.5*p4.y) + 
				(-0.5*p1.y + 0.5*p3.y);// + (p2.y);

	p.z = 3*pow(u,2)*(-0.5*p1.z + 1.5*p2.z - 1.5*p3.z + 0.5*p4.z) + 
				2*u*(p1.z - 2.5*p2.z + 2.0*p3.z - 0.5*p4.z) + 
				(-0.5*p1.z + 0.5*p3.z);// + (p2.z);

	p = unitize(p);
	return p;
}

point find_dderiv(double u, point p1, point p2, point p3, point p4)
{
	point p;
	
	p.x = 6*u*(-0.5*p1.x + 1.5*p2.x - 1.5*p3.x + 0.5*p4.x) + 
				2*(p1.x - 2.5*p2.x + 2.0*p3.x - 0.5*p4.x);
				// + (-0.5*p1.x + 0.5*p3.x);// + (p2.x);
	
	p.y = 6*u*(-0.5*p1.y + 1.5*p2.y - 1.5*p3.y + 0.5*p4.y) + 
				2*(p1.y - 2.5*p2.y + 2.0*p3.y - 0.5*p4.y);
				// + (-0.5*p1.y + 0.5*p3.y);// + (p2.y);

	p.z = 6*u*(-0.5*p1.z + 1.5*p2.z - 1.5*p3.z + 0.5*p4.z) + 
				2*(p1.z - 2.5*p2.z + 2.0*p3.z - 0.5*p4.z);
				// + (-0.5*p1.z + 0.5*p3.z);// + (p2.z);
	
	p = unitize(p);
	return p;
}

/* finds the tangent as 0.5*difference in the 
 * co-ordinates between points p1 and p2
 */
double find_dist(point p1, point p2)
{
	double dist;
	dist = sqrt(pow((p1.x - p2.x),2)+ pow((p1.y - p2.y),2) + 
							pow((p1.z - p2.z),2)); //REMEMBER: its square of the distance
	
	return dist;
}

void find_uRange(double u0, double u1, double maxLineLength, 
									point p1, point p2, point p3, point p4)
{
	double umid;
	double dist;

	umid = (u0+u1)/2;
	
	dist=find_dist(find_point(u0,p1,p2,p3,p4),find_point(u1,p1,p2,p3,p4));

	if(dist>maxLineLength)
	{
		find_uRange(u0,umid,maxLineLength,p1,p2,p3,p4);
		find_uRange(umid,u1,maxLineLength,p1,p2,p3,p4);
	}
	else
	{
		u.push_back(u0);
		u.push_back(u1);
	}
}

point crossProduct(point A,point B)
{
	point C;
	A=unitize(A);
	B=unitize(B);

	C.x = (A.y*B.z - B.y*A.z);
	C.y = (B.x*A.z - A.x*B.z);
	C.z = (A.x*B.y - A.y*B.x);
	
	C=unitize(C);
	return C;
}

void find_numOfSubSplines()
{
	for (int i=0;i<g_iNumOfSplines;i++)
	{
		g_iNumOfSubSplines+=g_Splines[i].numControlPoints;
	}
  g_subSplines = (struct spline *)malloc(g_iNumOfSubSplines * sizeof(struct spline));
  g_subTangents = (struct spline *)malloc(g_iNumOfSubSplines * sizeof(struct spline));
  g_subNormals = (struct spline *)malloc(g_iNumOfSubSplines * sizeof(struct spline));
  g_subBiNormals = (struct spline *)malloc(g_iNumOfSubSplines * sizeof(struct spline));
}

void generate_points()
{
	point p,q1,q2,q3,q4,startPoint,endPoint,exp;
	int numUValues;
	int subId;
	double derivMag;
	
	for (int i=0; i<g_iNumOfSplines; i++)
	{
		if(i==0)
		{
			startPoint.x = 0;
			startPoint.y = 0; 
			startPoint.z = 0;
		}
		else 
			startPoint=g_Splines[i-1].points[g_Splines[i-1].numControlPoints - 1];

		if(i==g_iNumOfSplines-1)
		{
			endPoint.x = 0;
			endPoint.y = 0; 
			endPoint.z = 0;
		}
		else	
			endPoint=g_Splines[i+1].points[0];
		
		exp.x=0;
		exp.y=1;
		exp.z=0;
	
		for (int j=0;j<g_Splines[i].numControlPoints-1;j++) 
		{
			q2=g_Splines[i].points[j];
			q3=g_Splines[i].points[j+1];
	
			if(j==0)
				q1 = startPoint;
			else	
				q1=g_Splines[i].points[j-1];
			
			if(j==g_Splines[i].numControlPoints-2)	
				q4=endPoint;
			else 
				q4=g_Splines[i].points[j+2];

			find_uRange(0,1,maxLineLength,q1,q2,q3,q4);
			numUValues=int(u.size());
			
			if (i>0)	
				subId=j+(g_Splines[i-1].numControlPoints);			
			else 
				subId = j;
			
			g_subSplines[subId].numControlPoints = numUValues;
			g_subTangents[subId].numControlPoints = numUValues;
			g_subNormals[subId].numControlPoints = numUValues;
			g_subBiNormals[subId].numControlPoints = numUValues;
			
			g_subSplines[subId].points= (struct	point *)malloc(numUValues * sizeof(struct point)); 
			g_subTangents[subId].points= (struct	point *)malloc(numUValues * sizeof(struct point)); 
			g_subNormals[subId].points= (struct	point *)malloc(numUValues * sizeof(struct point)); 
			g_subBiNormals[subId].points= (struct	point *)malloc(numUValues * sizeof(struct point)); 
			
 			g_subSplines[subId].points[0]=q2;
			g_subTangents[subId].points[0]=find_tangent(0,q1,q2,q3,q4);

			if(j==0)
				g_subNormals[subId].points[0] = exp;
			else
				g_subNormals[subId].points[0] = crossProduct(g_subBiNormals[subId-1].
																				points[g_subBiNormals[subId-1].numControlPoints-1], 
																				g_subTangents[subId].points[0]);
			
			g_subBiNormals[subId].points[0] = crossProduct(g_subTangents[subId].points[0], 
																				g_subNormals[subId].points[0]);

			for (int k=1;k<numUValues;k++)
			{
				g_subSplines[subId].points[k] = find_point(u[k-1],q1,q2,q3,q4); 
				g_subTangents[subId].points[k] = find_tangent(u[k-1],q1,q2,q3,q4);
				g_subNormals[subId].points[k] = crossProduct(g_subBiNormals[subId].points[k-1],
																				g_subTangents[subId].points[k]);
				g_subBiNormals[subId].points[k] = crossProduct(g_subTangents[subId].
																					points[k],g_subNormals[subId].points[k]);
			}
			u.clear();
		}
	}
}

void configTexture(int texId)
{
	// gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGB, my_sky.tWidth[texId], 
	// my_sky.tHeight[texId],GL_RGB,GL_UNSIGNED_BYTE,my_sky.sky_texels[texId]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,my_sky.tWidth[texId], 
			my_sky.tHeight[texId], 0, GL_RGB,GL_UNSIGNED_BYTE, 
			my_sky.sky_texels[texId]);
	
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
}

void drawSpline()
{
	int jhat,ihat;
	glBegin(GL_QUADS);
	glColor3f(0.5,0.5,0.3);
	for (int i=0;i<g_iNumOfSubSplines;i++)
	{
		jhat=20;
		ihat=i;
		for (int j=0;j<g_subSplines[i].numControlPoints-1;j+=20)
		{	
		//front
			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
								 	g_subSplines[i].points[j].y - 
								 	(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));
			
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x),
									g_subSplines[i].points[j].y - 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));
		//right	
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

		
		//left
			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y - 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y - 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

		
		//top	
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y - 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));
		
		//bottom
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y - 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

	
		//second spline
		
		//front
			glVertex3f(g_subSplines[i].points[j].x + 
									(dist*g_subBiNormals[i].points[j].x) - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(dist*g_subBiNormals[i].points[j].y) - 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z+(dist*g_subBiNormals[i].points[j].z) - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(dist*g_subBiNormals[i].points[j].x) + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(dist*g_subBiNormals[i].points[j].y) + 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z+(dist*g_subBiNormals[i].points[j].z) + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(dist*g_subBiNormals[i].points[j].x) + 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(dist*g_subBiNormals[i].points[j].y) + 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z+(dist*g_subBiNormals[i].points[j].z) + 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(dist*g_subBiNormals[i].points[j].x) - 
									(0.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(dist*g_subBiNormals[i].points[j].y) - 
									(0.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z+(dist*g_subBiNormals[i].points[j].z) - 
									(0.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

		//right	
			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

		//left
			glVertex3f(g_subSplines[i].points[j].x - 
									(-0.95*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y - 
									(-0.95*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z - 
									(-0.95*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(-0.95*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.95*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z+
									(0.95*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));  
			
			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

		//top	
			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.95*dist*g_subBiNormals[i].points[j].x) + 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.95*dist*g_subBiNormals[i].points[j].y) + 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(0.95*dist*g_subBiNormals[i].points[j].z) + 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) + 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].y) + 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) + 
									(0.0008*g_subNormals[ihat].points[jhat].z));

		//bottom
			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.95*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.95*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(0.95*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			
			if (j+20>=g_subSplines[i].numControlPoints-20)
			{
				if(i<g_iNumOfSubSplines-2) ihat=i+1;
				jhat=0;
			}
			else
				jhat=j+40;
		}
	}

	glEnd();

	glBegin(GL_QUADS);
	glColor3f(0.5,0.2,0.5);
	for (int i=0;i<g_iNumOfSubSplines;i+=4)
	{
		ihat=i;
		
		for (int j=70;j<g_subSplines[i].numControlPoints;j+=250)
		{
			jhat=j+2;
			//front
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y - 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));
			
			//back
			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			//left
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y),
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z + 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z + 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));
			
			//right
			glVertex3f(g_subSplines[i].points[j].x - 
								 (0.05*dist*g_subBiNormals[i].points[j].x) - 
								 (0.0008*g_subNormals[i].points[j].x), 
								 g_subSplines[i].points[j].y - 
								 (0.05*dist*g_subBiNormals[i].points[j].y) - 
								 (0.0008*g_subNormals[i].points[j].y), 
								 g_subSplines[i].points[j].z - 
								 (0.05*dist*g_subBiNormals[i].points[j].z) - 
								 (0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z - 
									(0.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[i].points[j].x - 
									(0.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z - 
									(0.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			//second spline	
			//front
			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.95*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.95*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.95*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(0.95*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z + 
									(0.95*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			//back	
			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			//left	
			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z + 
									(1.05*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.001*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[i].points[j].x + 
									(1.05*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(1.05*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(1.05*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));
									
		//right	
			glVertex3f(g_subSplines[i].points[j].x + 
									(0.95*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 
									g_subSplines[i].points[j].y + 
									(0.95*dist*g_subBiNormals[i].points[j].y) - 
									(0.0008*g_subNormals[i].points[j].y), 
									g_subSplines[i].points[j].z + 
									(0.95*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));
									
			glVertex3f(g_subSplines[i].points[j].x +
									(0.95*dist*g_subBiNormals[i].points[j].x) - 
									(0.0008*g_subNormals[i].points[j].x), 0.0, 
									g_subSplines[i].points[j].z + 
									(0.95*dist*g_subBiNormals[i].points[j].z) - 
									(0.0008*g_subNormals[i].points[j].z));
			
			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 0.0, 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

			glVertex3f(g_subSplines[ihat].points[jhat].x + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].x) - 
									(0.0008*g_subNormals[ihat].points[jhat].x), 
									g_subSplines[ihat].points[jhat].y + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].y) - 
									(0.0008*g_subNormals[ihat].points[jhat].y), 
									g_subSplines[ihat].points[jhat].z + 
									(0.95*dist*g_subBiNormals[ihat].points[jhat].z) - 
									(0.0008*g_subNormals[ihat].points[jhat].z));

		
		}
	}
	glEnd();
	
	glColor3f(1.0,1.0,1.0);	

	
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D,my_bar.crossTexture);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,my_bar.crossWidth,my_bar.crossHeight,0,GL_RGB,GL_UNSIGNED_BYTE,my_bar.cross_texels);
	
	glBegin(GL_QUADS);
	for (int i=0;i<g_iNumOfSubSplines;i++)
	{
		for (int j=0;j<g_subSplines[i].numControlPoints;j+=40)
		{
			glTexCoord2f(0, 1);
			glVertex3f(g_subSplines[i].points[j].x-(0.25*dist*g_subBiNormals[i].points[j].x),g_subSplines[i].points[j].y-(0.25*dist*g_subBiNormals[i].points[j].y),g_subSplines[i].points[j].z-(0.25*dist*g_subBiNormals[i].points[j].z));
      glTexCoord2f(1, 1); 
			glVertex3f(g_subSplines[i].points[j].x+(1.25*dist*g_subBiNormals[i].points[j].x),g_subSplines[i].points[j].y+(1.25*dist*g_subBiNormals[i].points[j].y),g_subSplines[i].points[j].z+(1.25*dist*g_subBiNormals[i].points[j].z));
      glTexCoord2f(1, 0);
			glVertex3f(g_subSplines[i].points[j+5].x+(1.25*dist*g_subBiNormals[i].points[j+5].x),g_subSplines[i].points[j+5].y+(1.25*dist*g_subBiNormals[i].points[j+5].y),g_subSplines[i].points[j+5].z+(1.25*dist*g_subBiNormals[i].points[j+5].z));
      glTexCoord2f(0, 0); 
			glVertex3f(g_subSplines[i].points[j+5].x-(0.25*dist*g_subBiNormals[i].points[j+5].x),g_subSplines[i].points[j+5].y-(0.25*dist*g_subBiNormals[i].points[j+5].y),g_subSplines[i].points[j+5].z-(0.25*dist*g_subBiNormals[i].points[j+5].z));
			
			glTexCoord2f(0, 1);
			glVertex3f(g_subSplines[i].points[j].x-(0.25*dist*g_subBiNormals[i].points[j].x),g_subSplines[i].points[j].y-(0.25*dist*g_subBiNormals[i].points[j].y),g_subSplines[i].points[j].z-(0.25*dist*g_subBiNormals[i].points[j].z));
      glTexCoord2f(1, 1); 
			glVertex3f(g_subSplines[i].points[j].x+(1.25*dist*g_subBiNormals[i].points[j].x),g_subSplines[i].points[j].y+(1.25*dist*g_subBiNormals[i].points[j].y),g_subSplines[i].points[j].z+(1.25*dist*g_subBiNormals[i].points[j].z));
      glTexCoord2f(1, 0);
			glVertex3f(g_subSplines[i].points[j].x+(1.25*dist*g_subBiNormals[i].points[j].x)-(0.1*dist*g_subNormals[i].points[j].x),g_subSplines[i].points[j].y+(1.25*dist*g_subBiNormals[i].points[j].y)-(0.1*dist*g_subNormals[i].points[j].y),g_subSplines[i].points[j].z+(1.25*dist*g_subBiNormals[i].points[j].z)-(0.1*dist*g_subNormals[i].points[j].z));
      glTexCoord2f(0, 0); 
			glVertex3f(g_subSplines[i].points[j].x-(0.25*dist*g_subBiNormals[i].points[j].x)-(0.1*dist*g_subNormals[i].points[j].x),g_subSplines[i].points[j].y-(0.25*dist*g_subBiNormals[i].points[j].y)-(0.1*dist*g_subNormals[i].points[j].y),g_subSplines[i].points[j].z-(0.25*dist*g_subBiNormals[i].points[j].z)-(0.1*dist*g_subNormals[i].points[j].z));
		}
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

void render_toy()
{
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D,my_toy.crossTexture);
		glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,my_toy.crossWidth,my_toy.crossHeight,0,GL_RGB,GL_UNSIGNED_BYTE,my_toy.cross_texels);

		glBegin(GL_QUADS);
		//front quad
		glTexCoord2f(0, 1); glVertex3f(  0.002f, -0.002f, -0.002f );
		glTexCoord2f(1, 1); glVertex3f( -0.002f, -0.002f, -0.002f );
		glTexCoord2f(1, 0); glVertex3f( -0.002f,  0.002f, -0.002f );
		glTexCoord2f(0, 0); glVertex3f(  0.002f,  0.002f, -0.002f );
		// Render the left quad
		glTexCoord2f(0, 1); glVertex3f(  0.002f, -0.002f,  0.002f );
		glTexCoord2f(1, 1); glVertex3f(  0.002f, -0.002f, -0.002f );
		glTexCoord2f(1, 0); glVertex3f(  0.002f,  0.002f, -0.002f );
		glTexCoord2f(0, 0); glVertex3f(  0.002f,  0.002f,  0.002f );
		// Render the right quad
		glTexCoord2f(0, 1); glVertex3f( -0.002f, -0.002f, -0.002f );
		glTexCoord2f(1, 1); glVertex3f( -0.002f, -0.002f,  0.002f );
		glTexCoord2f(1, 0); glVertex3f( -0.002f,  0.002f,  0.002f );
		glTexCoord2f(0, 0); glVertex3f( -0.002f,  0.002f, -0.002f );
		// Render the top quad
		glTexCoord2f(1, 1); glVertex3f( -0.002f,  0.002f, -0.002f );
		glTexCoord2f(1, 0); glVertex3f( -0.002f,  0.002f,  0.002f );
		glTexCoord2f(0, 0); glVertex3f(  0.002f,  0.002f,  0.002f );
		glTexCoord2f(0, 1); glVertex3f(  0.002f,  0.002f, -0.002f );
		// Render the bottom quad
		glTexCoord2f(0, 1); glVertex3f( -0.002f, -0.002f, -0.002f );
		glTexCoord2f(1, 1); glVertex3f( -0.002f, -0.002f,  0.002f );
		glTexCoord2f(1, 0); glVertex3f(  0.002f, -0.002f,  0.002f );
		glTexCoord2f(0, 0); glVertex3f(  0.002f, -0.002f, -0.002f );
		// Render the back quad
		glTexCoord2f(0, 1); glVertex3f( -0.002f, -0.002f,  0.002f );
		glTexCoord2f(1, 1); glVertex3f(  0.002f, -0.002f,  0.002f );
		glTexCoord2f(1, 0); glVertex3f(  0.002f,  0.002f,  0.002f );
		glTexCoord2f(0, 0); glVertex3f( -0.002f,  0.002f,  0.002f );
		glEnd();
		glDisable(GL_TEXTURE_2D);
}

void loadSky(int texId, char * filename)
{
	glGenTextures(1, &my_sky.skyTexture[texId]);

	texData = jpeg_read(filename, NULL);
	if (!texData)
	{
		printf ("error reading texture file.\n");
		exit(1);
	}
	
	my_sky.tWidth[texId] = texData->nx;
	my_sky.tHeight[texId] = texData->ny;
	int bpp = texData->bpp;
	my_sky.sky_texels[texId] = (GLubyte *)malloc(my_sky.tWidth[texId] * my_sky.tHeight[texId] * bpp * sizeof(GLubyte));
	memcpy(my_sky.sky_texels[texId],texData->pix,(my_sky.tWidth[texId] * my_sky.tHeight[texId] * bpp * sizeof(GLubyte)));
	//gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGB,my_sky.tWidth[texId],my_sky.tHeight[texId],GL_RGB,GL_UNSIGNED_BYTE,my_sky.sky_texels[texId]);
}

void loadTexture(int Id, char * filename)
{
	if (Id==0)
	{
		glGenTextures(1, &my_bar.crossTexture);
		texData = jpeg_read(filename, NULL);
		if (!texData)
		{
			printf ("error reading texture file.\n");
			exit(1);
		}

		my_bar.crossWidth = texData->nx;
		my_bar.crossHeight = texData->ny;
		int bpp = texData->bpp;
		my_bar.cross_texels = (GLubyte *)malloc(my_bar.crossWidth * my_bar.crossHeight* bpp * sizeof(GLubyte));
		memcpy(my_bar.cross_texels,texData->pix,(my_bar.crossWidth* my_bar.crossHeight* bpp * sizeof(GLubyte)));
	}
	else if(Id==1)
	{
		glGenTextures(1, &my_toy.crossTexture);
		texData = jpeg_read(filename, NULL);
		if (!texData)
		{
			printf ("error reading texture file.\n");
			exit(1);
		}
		
		my_toy.crossWidth = texData->nx;
		my_toy.crossHeight = texData->ny;
		int bpp = texData->bpp;
		my_toy.cross_texels = (GLubyte *)malloc(my_toy.crossWidth * my_toy.crossHeight* bpp * sizeof(GLubyte));
		memcpy(my_toy.cross_texels,texData->pix,(my_toy.crossWidth* my_toy.crossHeight* bpp * sizeof(GLubyte)));
	}
}


/* Textures */
void render_texture()
{
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_SRC_COLOR);

	// Front quad
  glBindTexture(GL_TEXTURE_2D,my_sky.skyTexture[0]);
  configTexture(0);
	glBegin(GL_QUADS);
      glTexCoord2f(0, 1); glVertex3f(  15.0f,		-5.0f, -15.0f );
      glTexCoord2f(1, 1); glVertex3f( -15.0f, 	-5.0f, -15.0f );
      glTexCoord2f(1, 0); glVertex3f( -15.0f,  15.0f, -15.0f );
      glTexCoord2f(0, 0); glVertex3f(  15.0f,  15.0f, -15.0f );
  glEnd();
	
	// Render the left quad
  glBindTexture(GL_TEXTURE_2D,my_sky.skyTexture[1]);
  configTexture(1);
  glBegin(GL_QUADS);
    glTexCoord2f(0, 1); glVertex3f(  15.0f, 	-5.0f,  15.0f );
    glTexCoord2f(1, 1); glVertex3f(  15.0f, 	-5.0f, -15.0f );
    glTexCoord2f(1, 0); glVertex3f(  15.0f,  15.0f, -15.0f );
		glTexCoord2f(0, 0); glVertex3f(  15.0f,  15.0f,  15.0f );
  glEnd();

	// Render the right quad
  glBindTexture(GL_TEXTURE_2D,my_sky.skyTexture[2]);
  configTexture(2);
  glBegin(GL_QUADS);
      glTexCoord2f(0, 1); glVertex3f( -15.0f, 	-5.0f, -15.0f );
      glTexCoord2f(1, 1); glVertex3f( -15.0f, 	-5.0f,  15.0f );
      glTexCoord2f(1, 0); glVertex3f( -15.0f,  15.0f,  15.0f );
      glTexCoord2f(0, 0); glVertex3f( -15.0f,  15.0f, -15.0f );
  glEnd();

	// Render the top quad
  glBindTexture(GL_TEXTURE_2D,my_sky.skyTexture[3]);
  configTexture(3);
  glBegin(GL_QUADS);
      glTexCoord2f(1, 1); glVertex3f( -15.0f,  15.0f, -15.0f );
      glTexCoord2f(1, 0); glVertex3f( -15.0f,  15.0f,  15.0f );
      glTexCoord2f(0, 0); glVertex3f(  15.0f,  15.0f,  15.0f );
      glTexCoord2f(0, 1); glVertex3f(  15.0f,  15.0f, -15.0f );
  glEnd();
 
	// Render the bottom quad
  glBindTexture(GL_TEXTURE_2D,my_sky.skyTexture[4]);
  configTexture(4);
  glBegin(GL_QUADS);
			glTexCoord2f(0, 1); glVertex3f( -15.0f, 	-5.0f, -15.0f );
      glTexCoord2f(1, 1); glVertex3f( -15.0f, 	-5.0f,  15.0f );
      glTexCoord2f(1, 0); glVertex3f(  15.0f, 	-5.0f,  15.0f );
      glTexCoord2f(0, 0); glVertex3f(  15.0f, 	-5.0f, -15.0f );
  glEnd();
	
	// Render the back quad
  glBindTexture(GL_TEXTURE_2D,my_sky.skyTexture[5]);
  configTexture(5);
  glBegin(GL_QUADS);
      glTexCoord2f(0, 1); glVertex3f( -15.0f, 	-5.0f,  15.0f );
      glTexCoord2f(1, 1); glVertex3f(  15.0f, 	-5.0f,  15.0f );
      glTexCoord2f(1, 0); glVertex3f(  15.0f,  15.0f,  15.0f );
      glTexCoord2f(0, 0); glVertex3f( -15.0f,  15.0f,  15.0f );
  glEnd();

	glDisable(GL_TEXTURE_2D);
	//glDeleteTextures(5, &my_sky.skyTexture[0]);
}

/* GL functions */
void myinit()
{
	glLineWidth(1.5f);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LINE_SMOOTH);
	glShadeModel(GL_SMOOTH);
}

/* create a window of size 640x480 */
void createWindow()
{
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH); 
	glutInitWindowSize (640,480); 
	glutInitWindowPosition (100,100);
	glutCreateWindow ("ShyamaDorbala_Assignment2");
	glEnable(GL_DEPTH_TEST);
}

/* reshape function to resize the window */
void reshape(int width, int height)
{
	glViewport(0,0,width,height); //setting the viewport to atleast contain the image
	glMatrixMode(GL_PROJECTION);// Projection Mode to view perspectively
  glLoadIdentity();// Initializing the Projection Matrix
	gluPerspective(60, (GLdouble) (width/height),0.01,1000); // Perspective Viewing with 60 degrees field of view
	
	glMatrixMode(GL_MODELVIEW); //getting back to model view
	glLoadIdentity(); // initializing it to start transformations
}

/*display callback */
void display()
{
	int bhat;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glLoadIdentity();
	
		glLoadIdentity();
		gluLookAt(g_subSplines[a].points[b].x+((dist/2)*g_subBiNormals[a].points[b].x)+0.01*g_subNormals[a].points[b].x,g_subSplines[a].points[b].y+((dist/2)*g_subBiNormals[a].points[b].y)+(0.01)*g_subNormals[a].points[b].y,g_subSplines[a].points[b].z+((dist/2)*g_subBiNormals[a].points[b].z)+(0.01)*g_subNormals[a].points[b].z,
		g_subSplines[a].points[b].x+((dist/2)*g_subBiNormals[a].points[b].x)+g_subTangents[a].points[b].x+0.02*g_subNormals[a].points[b].x,g_subSplines[a].points[b].y+((dist/2)*g_subBiNormals[a].points[b].y)+g_subTangents[a].points[b].y+0.02*g_subNormals[a].points[b].y,g_subSplines[a].points[b].z+((dist/2)*g_subBiNormals[a].points[b].z)+g_subTangents[a].points[b].z+0.02*g_subNormals[a].points[b].z,
		g_subNormals[a].points[b].x,g_subNormals[a].points[b].y,g_subNormals[a].points[b].z);
	/*else	
	{
		glTranslatef(0,-1,0);
		gluLookAt(0,0,-1,0,0,1,0,1,0); // Setting the camera position
	}*/

	glTranslatef(g_keyTranslate[0],g_keyTranslate[1],g_keyTranslate[2]);
	glTranslatef(g_vLandTranslate[0],g_vLandTranslate[1],g_vLandTranslate[2]);
	glRotatef(g_vLandRotate[0],1,0,0); //something is up with the axes of rotation
	glRotatef(g_vLandRotate[1],0,1,0);
	glRotatef(g_vLandRotate[2],0,0,1);
	glScalef(g_vLandScale[0],g_vLandScale[1],g_vLandScale[2]);
	
	render_texture();
	
	if(g_Sphere==1)
	{
		glPushMatrix();
		bhat=b+50;
		if(b<g_subSplines[a].numControlPoints-5)
			glTranslatef(g_subSplines[a].points[bhat].x+((dist/2)*g_subBiNormals[a].points[bhat].x),g_subSplines[a].points[bhat].y+((dist/2)*g_subBiNormals[a].points[bhat].y),g_subSplines[a].points[bhat].z+((dist/2)*g_subBiNormals[a].points[bhat].z));
		else if((a<g_iNumOfSubSplines-1) and ((bhat-g_subSplines[a].numControlPoints)<g_subSplines[a+1].numControlPoints))
			glTranslatef(g_subSplines[a+1].points[bhat-g_subSplines[a].numControlPoints].x+((dist/2)*g_subBiNormals[a+1].points[bhat-g_subSplines[a].numControlPoints].x),g_subSplines[a+1].points[bhat-g_subSplines[a].numControlPoints].y+((dist/2)*g_subBiNormals[a+1].points[bhat-g_subSplines[a].numControlPoints].y),g_subSplines[a+1].points[bhat-g_subSplines[a].numControlPoints].z+((dist/2)*g_subBiNormals[a+1].points[bhat-g_subSplines[a].numControlPoints].z));
		render_toy();
		glPopMatrix();
	}
	
	drawSpline();
	
	if(g_MoveCam==1)
	{
	if (a>=g_iNumOfSubSplines-2) a=0;
	else
	{
		if (g_subSplines[a+1].points[b].y-g_subSplines[a].points[b].y > 2)
			b+=1;
		else if (g_subSplines[a].points[b].y-g_subSplines[a+1].points[b].y > 0.1)
		{	
			b+=10; //FIXME: change this to control speed 
		}
		else
			b+=4;
	}
	}
//	if(b<g_subSplines[a].numControlPoints)
//	b+=4+round(sqrt(2*10*(maxHeight-g_subSplines[a].points[b].y)));
	//b+=4;
	if(b>=g_subSplines[a].numControlPoints-1) 
	{	
		b=0;
		a++;
		if (a>=g_iNumOfSubSplines-1) a=0;
	}
	glFlush();
	glutSwapBuffers();
}

void doIdle()
{
	string fname;
	static int k=0;
	
	glutPostRedisplay();
	// to take screenshot
	if (g_takeShot == 1) //comes from the menu
	{
		k++;

		if (k>30) 
		{
		stringstream ss;
    	ss << (k-30);
    	ss >> fname;
		fname+= ".jpg";
	
		if (k<40)
			fname.insert(0,"00");
		else if (k<130)
			fname.insert(0,"0");
	
		strcpy(fn, fname.c_str());

		 if (k<1030)
			 saveScreenshot(fn);
		//k+=1;
		}
	}
}

void mousebutton(int button, int state, int x, int y)
{
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
 
  switch(glutGetModifiers())
  {
    case GLUT_ACTIVE_CTRL:
      g_ControlState = TRANSLATE;
      break;
    case GLUT_ACTIVE_SHIFT:
      g_ControlState = SCALE;
      break;
    default:
      g_ControlState = ROTATE;
      break;
  }
}

void mouseidle(int x, int y)
{
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

/* converts mouse drags into information about 
rotation/translation/scaling */
void mousedrag(int x, int y)
{
  vMouseDelta[0] = x-g_vMousePos[0];
	vMouseDelta[1] = y-g_vMousePos[1];
  
  switch (g_ControlState)
  {
    case TRANSLATE:  
      if (g_iLeftMouseButton)
      {
        g_vLandTranslate[0] += vMouseDelta[0]*0.01;
        g_vLandTranslate[1] -= vMouseDelta[1]*0.01;
        g_vLandTranslate[2] += vMouseDelta[1]*0.01;
      }
      break;
    case ROTATE:
      if (g_iLeftMouseButton)
      {
        g_vLandRotate[0] += vMouseDelta[1];
        g_vLandRotate[1] += vMouseDelta[0];
        g_vLandRotate[2] += vMouseDelta[1];
      }
      break;
    case SCALE:
      if (g_iLeftMouseButton)
      {
        g_vLandScale[0] *= 1.0+vMouseDelta[0]*0.01;
        g_vLandScale[1] *= 1.0-vMouseDelta[1]*0.01;
        g_vLandScale[2] *= 1.0-vMouseDelta[1]*0.01;
      }
      break;
  }
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
} 

void keyboardFunc( unsigned char key, int x, int y )
{
	switch( key )
	{
    case 'q':
			exit( 1 );
    break;
		case 's':
			g_Sphere=1;
		break;
		case 't':
			g_takeShot=1;
		break;
		//case 'c':
		//	g_MoveCam=1;
		//	break;
  };

}
/* menu function for the gui */
void menufunc(int value)
{
	if (value == 3)
      exit(0);
	else if (value == 0)
			g_MoveCam = 0;
	else if (value == 1)
			g_MoveCam = 1;
	else if (value == 2)
			g_takeShot = 1;
	else if (value == 4)
	{
		g_Sphere=0;
		g_MoveCam=1;
	}
	else if (value == 5)
			g_Sphere=1;
}


void specialFunc( int key, int x, int y )
{
	int g_CtrlState=0;
	int g_ShiftState=0;

	if (glutGetModifiers()==GLUT_ACTIVE_CTRL)
      g_CtrlState = 1;
	else
			g_CtrlState = 0;
	
	if (glutGetModifiers()==GLUT_ACTIVE_SHIFT)
      g_ShiftState = 1;
	else
			g_ShiftState = 0;
	
	switch( key )
	{
		case GLUT_KEY_LEFT:
			g_keyTranslate[0]-=0.01;
			break;
		case GLUT_KEY_RIGHT:
			g_keyTranslate[0]+=0.01;
			break;
		case GLUT_KEY_UP:
			if (g_CtrlState==1)
				g_keyTranslate[2]+=0.01;
			else
				g_keyTranslate[1]+=0.01;
			break;
    case GLUT_KEY_DOWN:
			if (g_CtrlState==1)
				g_keyTranslate[2]-=0.01;
			else
				g_keyTranslate[1]-=0.01;
			break;
  };
}


int main (int argc, char ** argv)
{
  if (argc<2)
  {  
  printf ("usage: %s <trackfile> \n", argv[0]);
  exit(0);
  }

  loadSplines(argv[1]);
	find_numOfSubSplines();
	generate_points();
	maxHeight=find_maxHeight();
	
	loadSky(0,(char *) "./poszMed.jpg");
	loadSky(1,(char *) "./negxMed.jpg");
	loadSky(2,(char *) "./posxMed.jpg");
	loadSky(3,(char *) "./posyMed.jpg");
	loadSky(4,(char *) "./negyMed.jpg");
	loadSky(5,(char *) "./negzMed.jpg");
	loadTexture(0,(char *) "./bar.jpg"); 
	loadTexture(1,(char *) "./usc.jpg"); 

	glutInit(&argc,argv); 
	createWindow();
	myinit();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(doIdle);
	
	g_iMenuId = glutCreateMenu(menufunc);
  glutSetMenu(g_iMenuId);
  glutAddMenuEntry("Pause",0);
  glutAddMenuEntry("Play",1);
  glutAddMenuEntry("Add USC logo",5);
  glutAddMenuEntry("Remove USC logo",4);
  glutAddMenuEntry("Start taking screenshots",2);
  glutAddMenuEntry("Quit",3);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	
  glutMotionFunc(mousedrag);
  glutPassiveMotionFunc(mouseidle);
  glutMouseFunc(mousebutton);
	glutKeyboardFunc(keyboardFunc);
	glutSpecialFunc(specialFunc);

  glutMainLoop();
  return 0;

}
