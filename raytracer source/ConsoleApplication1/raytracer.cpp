#include "pch.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include "glm/glm.hpp"

#define MAX_DEPTH 3 // maximum depth for ray reflection

using namespace std;

void save_imageP6(int resx, int resy, char* fname, unsigned char* pixels);
void save_imageP3(int resx, int resy, char* fname, unsigned char* pixels);

// global variable (temporarily used for reading files)
float near;
float leftt;
float rightt;
float bottom;
float top;
float resx, resy;
string sname, lname, oname;
float stranx, strany, stranz, ssclx, sscly, ssclz, sr, sg, sb, ska, skd, sks, skr, sn;
float lposx, lposy, lposz, lr, lg, lb;
float br, bg, bb;
float ar, ag, ab;

// A class which store the properties of light sources
class Light
{
public:
	string lname;
	glm::vec3 lpos;
	glm::vec3 lrgb;

	// Constructor simply stores datas
	Light(string lnameT, glm::vec3 lposT, glm::vec3 lrgbT) : lname(lnameT), lpos(lposT), lrgb(lrgbT) {}
};

// A class which store the properties of sphere objects
class Sphere
{
public:
	string name;
	glm::vec3 tran; // tranlation part of transformation matrix
	glm::vec3 scl; // scaling part of transformation matrix
	float r;
	float g;
	float b;
	float ka;
	float kd;
	float ks;
	float kr;
	float n;

	// Following 2 matrixs are used to transform a ray in order to detect a intersection
	glm::mat4x4 objToWorld; // a matrix which transform transformation matrix of the sphere from object to world coordinate
	glm::mat4x4 worldToObj; // a matrix which transform transformation matrix of the sphere from world to object coordinate
	
	float radius;
	glm::vec3 sPos; // Sphere position, always (0,0,0) for this project because we are going to use "tran" values above to check intersection in world/object coordinates.

	// Constructor simply stores datas
	Sphere(
		string nameT,
		glm::vec3 tranT,
		glm::vec3 sclT,
		float rT,
		float gT,
		float bT,
		float kaT,
		float kdT,
		float ksT,
		float krT,
		float nT
	) : name(nameT),
		tran(tranT), scl(sclT), r(rT),
		g(gT), b(bT), ka(kaT), kd(kdT), ks(ksT), kr(krT), n(nT), radius(1.0f), sPos(glm::vec3(0.0f, 0.0f, 0.0f))
	{
		// an identity matrix to apply transformations
		glm::mat4x4 identity(1.0f);

		// translation needs to be multiplied prior to scaling
		glm::mat4x4 translation = {
		{ 1.0f, 0.0f, 0.0f, tranT.x },
		{ 0.0f, 1.0f, 0.0f, tranT.y },
		{ 0.0f, 0.0f, 1.0f, tranT.z },
		{ 0.0f, 0.0f, 0.0f, 1.0f }
		};

		objToWorld = translation * identity;

		glm::mat4x4 scale = {
		{ sclT.x, 0.0f, 0.0f, 0.0f },
		{ 0.0f, sclT.y, 0.0f, 0.0f },
		{ 0.0f, 0.0f, sclT.z, 0.0f },
		{ 0.0f, 0.0f, 0.0f, 1.0f }
		};

		objToWorld = scale * objToWorld;

		// transformation matrix that changes from world to object coordinate is simply inverse of objToWorld matrix
		worldToObj = glm::inverse(objToWorld);
	}

	// getter functions
	glm::mat4x4 getWorldToObj() { return worldToObj; }
	glm::mat4x4 getObjToWorld() { return objToWorld; }

	// determines if a ray intersects a sphere. 1st parameter is ray's starting point and 2nd parameter is direction of the ray.
	bool intersect(const glm::vec3 &rayOri, const glm::vec3 &rayDir, float &t0, float &t1) const
	{
		glm::vec3 center = rayOri - glm::vec3(0.0f, 0.0f, 0.0f); // sphere's center is always 0 as instructed in Assignment 3 pdf manual, but to move it later on the screen we use transformation matrix instead.
		float center_dot = glm::dot(center, center);
		float b = 2 * glm::dot(rayDir, center);

		float a = glm::dot(rayDir, rayDir);
		float c = center_dot - radius * radius;

		float root = (b * b) - (4 * a * c);
		if (root <= 0.0f) {
			return false;
		}

		t0 = (-b + sqrt(root)) / (2.0f * a);
		t1 = (-b - sqrt(root)) / (2.0f * a);

		// at least either one should be greater than near, so that we can render
		if (t0 >= near || t1 >= near) 
			return true;
		else
			return false;
	}

	// a function which checks light-sphere intersect
	bool intersectL(const glm::vec3 &rayOri, const glm::vec3 &rayDir, float &t0, float &t1)
	{
		glm::vec3 center = rayOri - glm::vec3(0.0f, 0.0f, 0.0f);
		float center_dot = glm::dot(center, center);
		float b = 2 * glm::dot(rayDir, center);

		float a = glm::dot(rayDir, rayDir);
		float c = center_dot - radius * radius;

		float root = b * b - 4 * a * c;
		if (root <= 0.f) {
			return false;
		}

		t0 = (-b + sqrt(root)) / (2.f * a);
		t1 = (-b - sqrt(root)) / (2.f * a);

		// at least either one should be greater than near so we can draw
		if (t0 >= near || t1 >= near) 
			return true;
		else
			return false;
	}
};

// 2 global variables that stores sphere and light objects
std::vector<Sphere> spheres;
std::vector<Light> lights;

// raytrace function
glm::vec3 raytrace(
	const glm::vec3 &rayOri, // ray origin in world coordinate
	const glm::vec3 &rayDir, // ray direction in world coordinate
	std::vector<Sphere> &spheres, // sphere objects
	const int &depth // for keeping depth count less than 3 times
)
{

	// return blank color if depth is greater than 3 for terminating reflection
	if (depth > MAX_DEPTH) {
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}

	glm::vec4 rayOriWld(rayOri.x, rayOri.y, rayOri.z, 1.0f); // a ray origin has w=1.0f because it is a point
	glm::vec4 rayDirWld(rayDir.x, rayDir.y, rayDir.z, 0.0f);// a ray direction has w=0.0f because it is a vector
	bool hollow = false;
	float tnear = 10000.0f;
	Sphere* sphere = NULL;
	glm::vec4 rayOriObj;
	glm::vec4 rayDirObj;

	// for each spheres
	for (unsigned i = 0; i < spheres.size(); ++i) {
		rayOriObj = rayOriWld * spheres[i].getWorldToObj(); // transform ray origin into object coordinate for checking intersection
		rayDirObj = rayDirWld * spheres[i].getWorldToObj(); // transform ray direction into object coordinate for checking intersection

		float t0 = 10000.0f, t1 = 10000.0f;
		if (spheres[i].intersect(rayOriObj, rayDirObj, t0, t1)) {
			
			// here, we store minimum possible because that is more close to our eye to render first
			if (t0 < t1 && t0 < tnear) {
				tnear = t0;
				sphere = &spheres[i];
			}
			else if (t1 < t0 && t1 < tnear) {
				tnear = t1;
				sphere = &spheres[i];
			}
		}

		// if intersecting sphere exists but it is hollow sphere, then we exchange the values of t0 and t1 because we need to render inside of the hollow sphere which corresponds to t1.
		if (sphere && tnear < near + 1 && depth == 0) {
			if (t0 > t1)
				tnear = t0;
			else
				tnear = t1;

			hollow = true;
		}

	}

	if (!sphere && depth > 0)
		return glm::vec3(0.0f, 0.0f, 0.0f); // if the reflected ray is not intersecting anything, return null color
	else if (!sphere)
		return glm::vec3(br, bg, bb); // return background color if the ray does not intersect anything

	///////////////////////////////////// COLORING - AMBIENT ///////////////////////////////////
	glm::vec3 totalColor(0.0f, 0.0f, 0.0f);
	glm::vec3 ambient(ar, ag, ab);
	glm::vec3 objColor(sphere->r, sphere->g, sphere->b);

	totalColor = sphere->ka * ambient * objColor;

	///////////////////////////////////// COLORING - PREPARATION ///////////////////////////////////

	// transform ray origin and direction to object space
	glm::vec4 rayOriObj2 = rayOriWld * sphere->getWorldToObj();
	glm::vec4 rayDirObj2 = rayDirWld * sphere->getWorldToObj();

	// calculate the point of intersection and its normal in object space
	glm::vec3 eyeToIntTimeObj = rayDirObj2 * tnear;
	glm::vec3 intersectionPointObj = glm::vec3(rayOriObj2.x, rayOriObj2.y, rayOriObj2.z) + eyeToIntTimeObj;
	glm::vec3 intersectingNormalObj = glm::normalize(intersectionPointObj - glm::vec3(0.0f, 0.0f, 0.0f)); // 오른쪽 이게 맞나 확인해보기

	// if the sphere is hollow, we need to negate the value because dot product of normals in inside of sphere and light direction must meet (which means, dot(A, B) < 0 so that they meet because their directions are not same)
	if (hollow)
		intersectingNormalObj *= -1;

	// transform intersection point, normal at intersection to world space
	glm::vec4 intersectionPointWld = glm::vec4(intersectionPointObj.x, intersectionPointObj.y, intersectionPointObj.z, 1.0f) * sphere->getObjToWorld();
	glm::mat4x4 wldToObj_transpose = glm::transpose(sphere->getWorldToObj());
	glm::vec4 intersectionNormalWld = glm::normalize(wldToObj_transpose * glm::vec4(intersectingNormalObj.x, intersectingNormalObj.y, intersectingNormalObj.z, 0.0f));
	
	// For each lights, we check if they intersect all objects for diffuse and specular calculation
	for (unsigned i = 0; i < lights.size(); ++i) {
		// if we have hollow sphere, deal with it first
		if (hollow) {
			// if hollow, we need to check if the light source is also in the inside of the sphere which only affects the inside of the hollow sphere
			float xSquare = lights[i].lpos.x * lights[i].lpos.x;
			float ySquare = lights[i].lpos.y * lights[i].lpos.y;
			float zSquare = lights[i].lpos.z * lights[i].lpos.z;
			float hollowRadius = sphere->radius * sphere->radius;

			// if satisfying below if statement, that means the light source is in inside of hollow sphere
			if ((xSquare + ySquare + zSquare) <= hollowRadius) {
				glm::vec3 dirFromInt2LightWld = glm::vec3(lights[i].lpos.x, lights[i].lpos.y, lights[i].lpos.z) - glm::vec3(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z);
				dirFromInt2LightWld = glm::normalize(dirFromInt2LightWld); // use it with intersectionPointWld to create a ray

				glm::vec4 intersectionPointObj;
				glm::vec4 dirFromInt2LightObj;
				float tnear2 = 10000.0f;
				Sphere* sphere2 = NULL;
				
				glm::vec3 dirFromLight2IntWld = -1.0f * dirFromInt2LightWld;
				if (glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), dirFromLight2IntWld) > 0) {
					// if the normal of intersection and direction of light has same direction, we do nothing for diffuse and specular because they do not meet ( "Do not meet" = "dot(A,B) > 0" = "have same direction" )
				}
				else {
					// For each light, check if there is intersection between light and sphere
					for (unsigned j = 0; j < spheres.size(); ++j) {

						intersectionPointObj = glm::vec4(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z, 1.0f) * spheres[j].getWorldToObj();
						dirFromInt2LightObj = glm::vec4(dirFromInt2LightWld.x, dirFromInt2LightWld.y, dirFromInt2LightWld.z, 0.0f) * spheres[j].getWorldToObj();
						float t0 = 10000.0f, t1 = 10000.0f;
					
						// Very similar to ray-sphere intersection part above
						if (spheres[j].intersectL(intersectionPointObj, dirFromInt2LightObj, t0, t1)) {
							if (t0 < t1 && t0 < tnear2) {
								tnear2 = t0;
								sphere2 = &spheres[j];
							}
							else if (t1 < t0 && t1 < tnear2) {
								tnear2 = t1;
								sphere2 = &spheres[j];
							}
						}

					}

					// If no intersection is found, apply local illumination. This also means the light is not blocked by others.
					if (!sphere2 || sphere == sphere2) { // if "sphere == sphere 2", that means intersection of eye-sphere and intersection of light-sphere happened in same sphere, which means, a hollow sphere in which the light source is inside of the sphere.
						// apply local illumination below

						///////////////////////////////////// COLORING - DIFFUSE ///////////////////////////////////
						
						// this part is already explained in our class slides and assignment 3 instruction pdf
						float NdotL = glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), dirFromInt2LightWld);
						NdotL = std::max(NdotL, 0.f);
						glm::vec3 secondEquation = sphere->kd * lights[i].lrgb * NdotL * glm::vec3(sphere->r, sphere->g, sphere->b);

						///////////////////////////////////// COLORING - SPECULAR ///////////////////////////////////
						
						// this part is already explained in our class slides and assignment 3 instruction pdf
						glm::vec3 dirFromInt2EyeWld = glm::normalize(glm::vec3(0.0f, 0.0f, 0.0f) - glm::vec3(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z)); // same as V
						float nDotL = glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), dirFromInt2LightWld);
						glm::vec3 R(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z);
						R *= 2 * NdotL;
						R -= dirFromInt2LightWld;
						float RdotV = glm::dot(R, dirFromInt2EyeWld);
						RdotV = std::max(RdotV, 0.f);

						glm::vec3 thirdEquation = sphere->ks * lights[i].lrgb * std::pow(RdotV, sphere->n);

						///////////////////////////////////// COLORING - TOTAL ADD ///////////////////////////////////
						totalColor += secondEquation;
						totalColor += thirdEquation;

					}
				} // end of for each light
			}
		} // end of hollow calculation
		else {
			// if not hollow sphere but regular sphere, deal with it. Calculation process is very similr to hollow sphere (I just separated them because just to prepare for case).
			// I am not going to explain this part because it is same as above which is already explained.
			glm::vec3 dirFromInt2LightWld = glm::vec3(lights[i].lpos.x, lights[i].lpos.y, lights[i].lpos.z) - glm::vec3(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z);
			dirFromInt2LightWld = glm::normalize(dirFromInt2LightWld); // use it with intersectionPointWld to create a ray

			glm::vec4 intersectionPointObj;
			glm::vec4 dirFromInt2LightObj;
			float tnear2 = 10000.0f;
			Sphere* sphere2 = NULL;
			
			glm::vec3 dirFromLight2IntWld = -1.0f * dirFromInt2LightWld;
			if (glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), dirFromLight2IntWld) > 0) {
			}
			else {
				for (unsigned j = 0; j < spheres.size(); ++j) {

					intersectionPointObj = glm::vec4(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z, 1.0f) * spheres[j].getWorldToObj();
					dirFromInt2LightObj = glm::vec4(dirFromInt2LightWld.x, dirFromInt2LightWld.y, dirFromInt2LightWld.z, 0.0f) * spheres[j].getWorldToObj();
					float t0 = 10000.0f, t1 = 10000.0f;

					if (spheres[j].intersectL(intersectionPointObj, dirFromInt2LightObj, t0, t1)) {
						if (t0 < t1 && t0 < tnear2) {
							tnear2 = t0;
							sphere2 = &spheres[j];
						}
						else if (t1 < t0 && t1 < tnear2) {
							tnear2 = t1;
							sphere2 = &spheres[j];
						}
					}

				}

				// This part is also same as above except that it is not for light-hollow sphere but light-regular sphere intersection. 
				// I am not going to explain this part because it is same as above which is already explained.
				if (!sphere2) {
					///////////////////////////////////// COLORING - DIFFUSE ///////////////////////////////////
					float NdotL = glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), dirFromInt2LightWld);
					NdotL = std::max(NdotL, 0.f);

					glm::vec3 secondEquation = sphere->kd * lights[i].lrgb * NdotL * glm::vec3(sphere->r, sphere->g, sphere->b);

					///////////////////////////////////// COLORING - SPECULAR ///////////////////////////////////
					glm::vec3 dirFromInt2EyeWld = glm::normalize(glm::vec3(0.0f, 0.0f, 0.0f) - glm::vec3(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z)); // same as V
					float nDotL = glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), dirFromInt2LightWld);
					glm::vec3 R(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z);
					R *= 2 * NdotL;
					R -= dirFromInt2LightWld;
					float RdotV = glm::dot(R, dirFromInt2EyeWld);
					RdotV = std::max(RdotV, 0.f);

					glm::vec3 thirdEquation = sphere->ks * lights[i].lrgb * std::pow(RdotV, sphere->n);

					///////////////////////////////////// COLORING - TOTAL ADD ///////////////////////////////////
					totalColor += secondEquation;
					totalColor += thirdEquation;

				}
			}
		}
	} // end of for each light

	///////////////////////////////////// COLORING - REFLECTION ///////////////////////////////////
	if (sphere->kr > 0.0f) {
		// How to get reflected ray is explained in our lecture slide and following code does same as explained there
		glm::vec3 c(rayDirWld.x, rayDirWld.y, rayDirWld.z);
		glm::vec3 p(intersectionPointWld.x, intersectionPointWld.y, intersectionPointWld.z);
		float Ndotc = glm::dot(glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z), c);
		glm::vec3 v = (-2 * Ndotc * glm::vec3(intersectionNormalWld.x, intersectionNormalWld.y, intersectionNormalWld.z)) + c;

		glm::vec3 tempColor(0.0f, 0.0f, 0.0f);
		tempColor = raytrace(p, v, spheres, depth + 1); // retrace with new ray for reflection, while giving depth + 1 value to prevent reflecting more than 3 times
		tempColor *= sphere->kr;
		totalColor += tempColor;
	}

	// return our final color
	return totalColor;

}


void render(std::vector<Sphere> &spheres)
{

	glm::vec3 pixel; // a pixel to store color returnbed by raytrace function

	char fname3[20] = "sceneP3.ppm"; //This should be set based on the input file
	char fname6[20] = "sceneP6.ppm"; //This should be set based on the input file
	unsigned char *pixels;

	// This will be your image. Note that pixels[0] is the top left of the image and
	// pixels[3*resx*resy-1] is the bottom right of the image.
	pixels = new unsigned char[3 * resx*resy];

	int k = 0;

	// for each col and row of the screen
	for (unsigned y = 0; y < resy; y++) {
		for (unsigned x = 0; x < resx; x++) {
			// a ray direction in world coordinates. if x=0 and y=0, then current pixel will correspond to (left, -botom) in world coordinate rendered at bottom-left conner of the image.
			// similarly, if x=width, y=height, then current pixel will correspond to (right, top) in world coordinate rendered at top-right conner of the image.
			float xx = leftt + (rightt * 2.0f * (float)x / (float)resx);
			float yy = -bottom + (-top * 2.0f * (float)y / (float)resy);

			glm::vec3 rayDir(xx, yy, -1);
			// rayDir = glm::normalize(rayDir);

			// generate ray which has origin(eye) at (0,0,0) and direction to world coordinate of the pixel, and trace.
			pixel = raytrace(glm::vec3(0.0f, 0.0f, 0.0f), rayDir, spheres, 0);

			// sometimes values can go higher than its maximum. Fix to 1.0 if it happens.
			if (pixel.x > 1.0f)
				pixel.x = 1.0f;
			if (pixel.y > 1.0f)
				pixel.y = 1.0f;
			if (pixel.z > 1.0f)
				pixel.z = 1.0f;

			// to render in PPM file
			pixels[k] = pixel.x * 255;
			pixels[k + 1] = pixel.y * 255;
			pixels[k + 2] = pixel.z * 255;

			k = k + 3;
		}
	}

	char *cstr = &oname[0u];
	//save_imageP3(resx, resy, cstr, pixels);
	save_imageP6(resx, resy, cstr, pixels); // Your final program must produce the binary version (P6).

}


bool readValues(std::stringstream &s, const int numValues, float* values)
{
	for (int i = 0; i < numValues; i++) {
		s >> values[i]; // store values into temporary float array
		if (s.fail()) { // if any errors, return false
			return false;
		}
	}
	return true; // return true if all datas from the file in a line is stored in temporary float array
}

// a function that deal with a file
void readFile(const char *fileName) {

	string line, first;
	ifstream in;
	in.open(fileName); // open a file
	if (in.is_open()) {
		getline(in, line); // read a line
		while (in) {
			stringstream stream(line);
			stream >> first; // read first word in a line
			float values[14];

			bool isValid;

			// if first word is NEAR
			if (first == "NEAR") {
				isValid = readValues(stream, 1, values); // temporarily store its values in "float values[14]"
				if (isValid) { // if successfully read all datas in a line and they are stored in "float values[14]", we store in our actual variables.
					near = values[0];
				}
			}

			else if (first == "LEFT") {
				isValid = readValues(stream, 1, values);
				if (isValid) {
					leftt = values[0];
				}
			}

			else if (first == "RIGHT") {
				isValid = readValues(stream, 1, values);
				if (isValid) {
					rightt = values[0];
				}
			}

			else if (first == "BOTTOM") {
				isValid = readValues(stream, 1, values);
				if (isValid) {
					bottom = values[0];
				}
			}

			else if (first == "TOP") {
				isValid = readValues(stream, 1, values);
				if (isValid) {
					top = values[0];
				}
			}

			else if (first == "OUTPUT") {
				stream >> oname;
			}

			else if (first == "RES") {
				isValid = readValues(stream, 2, values);
				if (isValid) {
					resx = (float)values[0];
					resy = (float)values[1];
				}
			}

			else if (first == "BACK") {
				isValid = readValues(stream, 3, values);
				if (isValid) {
					br = (float)values[0];
					bg = (float)values[1];
					bb = (float)values[2];
				}
			}

			else if (first == "AMBIENT") {
				isValid = readValues(stream, 3, values);
				if (isValid) {
					ar = (float)values[0];
					ag = (float)values[1];
					ab = (float)values[2];
				}
			}

			else if (first == "SPHERE") {
				stream >> sname;
				isValid = readValues(stream, 14, values);
				if (isValid) {
					stranx = (float)values[0];
					strany = (float)values[1];
					stranz = (float)values[2];
					ssclx = (float)values[3];
					sscly = (float)values[4];
					ssclz = (float)values[5];
					sr = (float)values[6];
					sg = (float)values[7];
					sb = (float)values[8];
					ska = (float)values[9];
					skd = (float)values[10];
					sks = (float)values[11];
					skr = (float)values[12];
					sn = (float)values[13];

					spheres.push_back(Sphere(sname, glm::vec3(stranx, strany, stranz), glm::vec3(ssclx, sscly, ssclz), sr, sg, sb, ska, skd, sks, skr, sn));
				}
			}

			else if (first == "LIGHT") {
				stream >> lname;
				isValid = readValues(stream, 6, values);
				if (isValid) {
					lposx = (float)values[0];
					lposy = (float)values[1];
					lposz = (float)values[2];
					lr = (float)values[3];
					lg = (float)values[4];
					lb = (float)values[5];

					lights.push_back(Light(lname, glm::vec3(lposx, lposy, lposz), glm::vec3(lr, lg, lb)));
				}
			}
			getline(in, line);
		}

	}
	in.close(); // close the file
}



int main(int argc, char **argv)
{
	// we read file to store values there
	readFile(argv[1]);

	// then, we render
	render(spheres);

	return 0;
}

// Output in P6 format, a binary file containing:
// P6
// ncolumns nrows
// Max colour value
// colours in binary format thus unreadable
void save_imageP6(int resx, int resy, char* fname, unsigned char* pixels) {
	FILE *fp;
	const int maxVal = 255;

	printf("Saving image %s: %d x %d\n", fname, resx, resy);
	//fp = fopen_s(fname, "wb");
	fopen_s(&fp, fname, "wb");
	if (!fp) {
		printf("Unable to open file '%s'\n", fname);
		return;
	}
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", resx, resy);
	fprintf(fp, "%d\n", maxVal);

	for (int j = 0; j < resy; j++) {
		fwrite(&pixels[j*resx * 3], 3, resx, fp);
	}

	fclose(fp);
}

// Output in P3 format, a text file containing:
// P3
// ncolumns nrows
// Max colour value (for us, and usually 255)
// r1 g1 b1 r2 g2 b2 .....
void save_imageP3(int resx, int resy, char* fname, unsigned char* pixels) {
	FILE *fp;
	const int maxVal = 255;

	printf("Saving image %s: %d x %d\n", fname, resx, resy);
	//fp = fopen(fname, "w");
	fopen_s(&fp, fname, "w");
	if (!fp) {
		printf("Unable to open file '%s'\n", fname);
		return;
	}
	fprintf(fp, "P3\n");
	fprintf(fp, "%d %d\n", resx, resy);
	fprintf(fp, "%d\n", maxVal);

	int k = 0;
	for (int j = 0; j < resy; j++) {

		for (int i = 0; i < resx; i++)
		{
			fprintf(fp, " %d %d %d", pixels[k], pixels[k + 1], pixels[k + 2]);
			k = k + 3;
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
