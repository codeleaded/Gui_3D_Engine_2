#include <math.h>
#include "C:/Wichtig/System/Static/Library/WindowEngine1.0.h"
#include "C:\Wichtig\System\Static\Container\Vector.h"

typedef struct vec3d
{
	float x, y, z;
} vec3d;

typedef struct triangle
{
	vec3d p[3];
	unsigned int col;
} triangle;

typedef struct mesh
{
	Vector tris;
} mesh;

typedef struct mat4x4
{
	float m[4][4];
} mat4x4;

mesh meshCube;
mat4x4 matProj;

vec3d vCamera;

double fTheta;

void MultiplyMatrixVector(vec3d i, vec3d* o, mat4x4 m){
	o->x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
	o->y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
	o->z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
	float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];
	if (w != 0.0f)
	{
		o->x /= w; o->y /= w; o->z /= w;
	}
}

int QCompare(const void* p1,const void* p2){
	triangle* t1 = (triangle*)p1;
	triangle* t2 = (triangle*)p2;
	float z1 = (t1->p[0].z + t1->p[1].z + t1->p[2].z) / 3.0f;
	float z2 = (t2->p[0].z + t2->p[1].z + t2->p[2].z) / 3.0f;
	return z1 > z2;
}

void Setup(AlxWindow* w){
    meshCube.tris = Vector_New(sizeof(triangle));
    
    // SOUTH                                                     
    triangle tri = {0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f};
    Vector_Push(&meshCube.tris,&tri);
    tri = (triangle){0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f};
    Vector_Push(&meshCube.tris,&tri);
    
    //EAST                                                     
    tri = (triangle){ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (triangle){ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // NORTH                                                     
    tri = (triangle){ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (triangle){ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // WEST                                                      
    tri = (triangle){ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (triangle){ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // TOP                                                       
    tri = (triangle){ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (triangle){ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    
    // BOTTOM                                                    
    tri = (triangle){ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);
    tri = (triangle){ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f };
    Vector_Push(&meshCube.tris,&tri);


	float fNear = 0.1f;
	float fFar = 1000.0f;
	float fFov = 90.0f;
	float fAspectRatio = (float)GetHeight() / (float)GetWidth();
	float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);
	memset(matProj.m[0],0,4 * sizeof(float));
    memset(matProj.m[1],0,4 * sizeof(float));
    memset(matProj.m[2],0,4 * sizeof(float));
    memset(matProj.m[3],0,4 * sizeof(float));
    matProj.m[0][0] = fAspectRatio * fFovRad;
	matProj.m[1][1] = fFovRad;
	matProj.m[2][2] = fFar / (fFar - fNear);
	matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matProj.m[2][3] = 1.0f;
	matProj.m[3][3] = 0.0f;
}

void Update(AlxWindow* w){
    // Clear Screen
	Clear(0);
	// Set up rotation matrices
	mat4x4 matRotZ, matRotX;
	fTheta += 1.0 * w->ElapsedTime;
	
    // Rotation Z
    memset(matRotZ.m[0],0,4 * sizeof(float));
    memset(matRotZ.m[1],0,4 * sizeof(float));
    memset(matRotZ.m[2],0,4 * sizeof(float));
    memset(matRotZ.m[3],0,4 * sizeof(float));
	matRotZ.m[0][0] = cosf(fTheta);
	matRotZ.m[0][1] = sinf(fTheta);
	matRotZ.m[1][0] = -sinf(fTheta);
	matRotZ.m[1][1] = cosf(fTheta);
	matRotZ.m[2][2] = 1;
	matRotZ.m[3][3] = 1;
	// Rotation X
    memset(matRotX.m[0],0,4 * sizeof(float));
    memset(matRotX.m[1],0,4 * sizeof(float));
    memset(matRotX.m[2],0,4 * sizeof(float));
    memset(matRotX.m[3],0,4 * sizeof(float));
	matRotX.m[0][0] = 1;
	matRotX.m[1][1] = cosf(fTheta * 0.5f);
	matRotX.m[1][2] = sinf(fTheta * 0.5f);
	matRotX.m[2][1] = -sinf(fTheta * 0.5f);
	matRotX.m[2][2] = cosf(fTheta * 0.5f);
	matRotX.m[3][3] = 1;

	Vector vecTrianglesToRaster = Vector_New(sizeof(triangle));
	
	for (int i = 0;i<meshCube.tris.size;i++){
        triangle tri = *(triangle*)Vector_Get(&meshCube.tris,i);
		triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;
		// Rotate in Z-Axis
		MultiplyMatrixVector(tri.p[0], &triRotatedZ.p[0], matRotZ);
		MultiplyMatrixVector(tri.p[1], &triRotatedZ.p[1], matRotZ);
		MultiplyMatrixVector(tri.p[2], &triRotatedZ.p[2], matRotZ);
		// Rotate in X-Axis
		MultiplyMatrixVector(triRotatedZ.p[0], &triRotatedZX.p[0], matRotX);
		MultiplyMatrixVector(triRotatedZ.p[1], &triRotatedZX.p[1], matRotX);
		MultiplyMatrixVector(triRotatedZ.p[2], &triRotatedZX.p[2], matRotX);
		// Offset into the screen
		triTranslated = triRotatedZX;
		triTranslated.p[0].z = triRotatedZX.p[0].z + 3.0f;
		triTranslated.p[1].z = triRotatedZX.p[1].z + 3.0f;
		triTranslated.p[2].z = triRotatedZX.p[2].z + 3.0f;

		vec3d normal, line1, line2;
		line1.x = triTranslated.p[1].x - triTranslated.p[0].x;
		line1.y = triTranslated.p[1].y - triTranslated.p[0].y;
		line1.z = triTranslated.p[1].z - triTranslated.p[0].z;

		line2.x = triTranslated.p[2].x - triTranslated.p[0].x;
		line2.y = triTranslated.p[2].y - triTranslated.p[0].y;
		line2.z = triTranslated.p[2].z - triTranslated.p[0].z;

		normal.x = line1.y * line2.z - line1.z * line2.y;
		normal.y = line1.z * line2.x - line1.x * line2.z;
		normal.z = line1.x * line2.y - line1.y * line2.x;
		float l = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
		normal.x /= l; normal.y /= l; normal.z /= l;

		//if (normal.z < 0)
		if(normal.x * (triTranslated.p[0].x - vCamera.x) + 
		   normal.y * (triTranslated.p[0].y - vCamera.y) +
		   normal.z * (triTranslated.p[0].z - vCamera.z) < 0.0f)
		{
			vec3d light_direction = { 0.0f, 0.0f, -1.0f };
			float l = sqrtf(light_direction.x*light_direction.x + light_direction.y*light_direction.y + light_direction.z*light_direction.z);
			light_direction.x /= l; light_direction.y /= l; light_direction.z /= l;
			
			float dp = normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z;
			
			triTranslated.col = dp * 255;

			MultiplyMatrixVector(triTranslated.p[0], &triProjected.p[0], matProj);
			MultiplyMatrixVector(triTranslated.p[1], &triProjected.p[1], matProj);
			MultiplyMatrixVector(triTranslated.p[2], &triProjected.p[2], matProj);
			triProjected.col = triTranslated.col;
			
			triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
			triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
			triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;
			triProjected.p[0].x *= 0.5f * (float)GetWidth();
			triProjected.p[0].y *= 0.5f * (float)GetHeight();
			triProjected.p[1].x *= 0.5f * (float)GetWidth();
			triProjected.p[1].y *= 0.5f * (float)GetHeight();
			triProjected.p[2].x *= 0.5f * (float)GetWidth();
			triProjected.p[2].y *= 0.5f * (float)GetHeight();
			
			Vector_Push(&vecTrianglesToRaster,&triProjected);
		}
	}

	qsort(vecTrianglesToRaster.Memory,vecTrianglesToRaster.size,vecTrianglesToRaster.ELEMENT_SIZE,QCompare);

	for(int i = 0;i<vecTrianglesToRaster.size;i++){
		triangle triProjected = *(triangle*)Vector_Get(&vecTrianglesToRaster,i);
		//RenderTriangle((Vec2){triProjected.p[0].x,triProjected.p[0].y},
		//	                   (Vec2){triProjected.p[1].x,triProjected.p[1].y},
		//	                   (Vec2){triProjected.p[2].x,triProjected.p[2].y},
		//	                   (Pixel){triProjected.sym,triProjected.col});
		RenderTriangleWire(((Vec2){triProjected.p[0].x,triProjected.p[0].y}),
			               ((Vec2){triProjected.p[1].x,triProjected.p[1].y}),
			               ((Vec2){triProjected.p[2].x,triProjected.p[2].y}),
			               (Pixel){0xFFFFFFFF},1.0f);
	}
	Vector_Free(&vecTrianglesToRaster);
}

void Delete(AlxWindow* w){
    Vector_Free(&meshCube.tris);
}

int main(){
    if(Create("Game Test",400,300,4,4,Setup,Update,Delete))
        Start();
    return 0;
}