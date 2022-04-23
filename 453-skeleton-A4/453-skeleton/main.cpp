#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"
#include "Camera.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

#define _USE_MATH_DEFINES
#include <math.h>




// An example struct for Game Objects.
// You are encouraged to customize this as you see fit.
struct GameObject {
	// Struct's constructor deals with the texture.
	// Also sets default position, theta, scale, and transformationMatrix
	GameObject(std::string texturePath, GLenum textureInterpolation) :
		texture(texturePath, textureInterpolation),
		position(0.0f, 0.0f, 0.0f),
		theta(glm::radians(90.f)),
		scale(1),
		transformationMatrix(1.0f) // This constructor sets it as the identity matrix
	{}

	CPU_Geometry cgeom;
	GPU_Geometry ggeom;
	Texture texture;

	glm::vec3 position;
	float theta; // Object's rotation
	// Alternatively, you could represent rotation via a normalized heading vec:
	// glm::vec3 heading;
	glm::vec3 scale; // Or, alternatively, a glm::vec2 scale;
	glm::mat4 transformationMatrix;
};

// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
	gpuGeom.setNormals(cpuGeom.normals);
	gpuGeom.setTexCoords(cpuGeom.texCoords);
}

// EXAMPLE CALLBACKS
class Assignment4 : public CallbackInterface {

public:
	Assignment4() :
		camera(0.0, 0.0, 2.0),
		aspect(1.0f),
		spacePressed(false),
		rPressed(false)
	{}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
			spacePressed = !spacePressed;
		}
		//if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE) {
		//	spacePressed = false;
		//}
		if (key == GLFW_KEY_R && action == GLFW_PRESS) {
			rPressed = true;
		}
		if (key == GLFW_KEY_R && action == GLFW_RELEASE) {
			rPressed = false;
		}
	}

	bool buttonDown(int button) {
		if (button == GLFW_KEY_SPACE) {
			return spacePressed;
		}
		if (button == GLFW_KEY_R) {
			return rPressed;
		}
		return false;
	}

	void refreshStatuses() {
		//spacePressed = false;
		rPressed = false;
	}

	virtual void mouseButtonCallback(int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_RIGHT) {
			if (action == GLFW_PRESS) {
				rightMouseDown = true;
			} else if (action == GLFW_RELEASE) {
				rightMouseDown = false;
			}
		}
	}
	virtual void cursorPosCallback(double xpos, double ypos) {
		if (rightMouseDown) {
			double dx = xpos - mouseOldX;
			double dy = ypos - mouseOldY;
			camera.incrementTheta(dy);
			camera.incrementPhi(dx);
		}
		mouseOldX = xpos;
		mouseOldY = ypos;
	}
	virtual void scrollCallback(double xoffset, double yoffset) {
		camera.incrementR(yoffset);
	}
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
		aspect = float(width)/float(height);
	}

	void viewPipeline(ShaderProgram &sp) {
		glm::mat4 M = glm::mat4(1.0);
		glm::mat4 V = camera.getView();
		glm::mat4 P = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 1000.f);

		GLint location = glGetUniformLocation(sp, "light");
		glm::vec3 light = camera.getPos();
		glUniform3fv(location, 1, glm::value_ptr(light));

		GLint uniMat = glGetUniformLocation(sp, "M");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(M));
		uniMat = glGetUniformLocation(sp, "V");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(V));
		uniMat = glGetUniformLocation(sp, "P");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(P));
	}



	Camera camera;

private:

	bool rightMouseDown = false;
	float aspect;
	double mouseOldX;
	double mouseOldY;
	bool spacePressed = false;
	bool rPressed = false;

};

void colouredTriangles(CPU_Geometry &geom) {
	geom.cols.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.cols.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.cols.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.cols.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.cols.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.cols.push_back(glm::vec3(1.0, 0.0, 0.0));
}

void positiveZFace(std::vector<glm::vec3> const &originQuad, CPU_Geometry &geom) {
	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, 0.5));
	for(auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(
			glm::vec3(T * glm::vec4((*i), 1.0))
		);
	}
	geom.normals.push_back(glm::vec3(0.0, 0.0, 1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, 1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, 1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, 1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, 1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, 1.0));
}

void positiveXFace(std::vector<glm::vec3> const &originQuad, CPU_Geometry &geom) {
	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, 0.0f, 0.0f));
	for(auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(
			glm::vec3(T * R * glm::vec4((*i), 1.0))
		);
	}
	geom.normals.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(1.0, 0.0, 0.0));
}

void negativeZFace(std::vector<glm::vec3> const &originQuad, CPU_Geometry &geom) {
	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -0.5f));
	for(auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(
			glm::vec3(T * R * glm::vec4((*i), 1.0))
		);
	}
	geom.normals.push_back(glm::vec3(0.0, 0.0, -1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, -1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, -1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, -1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, -1.0));
	geom.normals.push_back(glm::vec3(0.0, 0.0, -1.0));
}

void negativeXFace(std::vector<glm::vec3> const &originQuad, CPU_Geometry &geom) {
	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, 0.0f, 0.0f));
	for(auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(
			glm::vec3(T * R * glm::vec4((*i), 1.0))
		);
	}
	geom.normals.push_back(glm::vec3(-1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(-1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(-1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(-1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(-1.0, 0.0, 0.0));
	geom.normals.push_back(glm::vec3(-1.0, 0.0, 0.0));
}

void positiveYFace(std::vector<glm::vec3> const &originQuad, CPU_Geometry &geom) {
	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.5f, 0.0f));
	for(auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(
			glm::vec3(T * R * glm::vec4((*i), 1.0))
		);
	}
	geom.normals.push_back(glm::vec3(0.0, 1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, 1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, 1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, 1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, 1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, 1.0, 0.0));
}

void negativeYFace(std::vector<glm::vec3> const &originQuad, CPU_Geometry &geom) {
	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -0.5f, 0.0f));
	for(auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(
			glm::vec3(T * R * glm::vec4((*i), 1.0))
		);
	}
	geom.normals.push_back(glm::vec3(0.0, -1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, -1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, -1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, -1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, -1.0, 0.0));
	geom.normals.push_back(glm::vec3(0.0, -1.0, 0.0));
}

void createSphere(float radius, int degree, CPU_Geometry& sphere) {
	float increment = M_PI / degree;
	radius = 0.5;

	int x = 1;

	for (float phi = 0; phi < 1 * M_PI; phi += increment) {
		for (float theta = 0; theta < 2 * M_PI; theta += increment) {
			sphere.verts.push_back(glm::vec3{ radius * cos(theta) * sin(phi), radius * sin(theta) * sin(phi), radius * cos(phi) });
			sphere.normals.push_back(glm::vec3{ cos(theta) * sin(phi),  sin(theta) * sin(phi), cos(phi) });
			sphere.texCoords.push_back(glm::vec2{1- theta / (2 * M_PI), phi / (1 * M_PI) });

			sphere.verts.push_back(glm::vec3{ radius * cos(theta + increment) * sin(phi), radius * sin(theta + increment) * sin(phi), radius * cos(phi) });
			sphere.normals.push_back(glm::vec3{ cos(theta + increment) * sin(phi),  sin(theta + increment) * sin(phi), cos(phi) });
			sphere.texCoords.push_back (glm::vec2{1- (theta + increment) / (2 * M_PI), phi / (1 * M_PI) });

			sphere.verts.push_back(glm::vec3{ radius * cos(theta) * sin(phi + increment), radius * sin(theta) * sin(phi + increment), radius * cos(phi + increment) });
			sphere.normals.push_back(glm::vec3{ cos(theta) * sin(phi + increment),  sin(theta) * sin(phi + increment),  cos(phi + increment) });
			sphere.texCoords.push_back(glm::vec2{1- theta / (2 * M_PI), (phi + increment) / (1 * M_PI) });

			sphere.verts.push_back(glm::vec3{ radius * cos(theta) * sin(phi + increment), radius * sin(theta) * sin(phi + increment), radius * cos(phi + increment) });
			sphere.normals.push_back(glm::vec3{ cos(theta) * sin(phi + increment), sin(theta) * sin(phi + increment), cos(phi + increment) });
			sphere.texCoords.push_back(glm::vec2{1- theta / (2 * M_PI), (phi + increment) / (1 * M_PI) });

			sphere.verts.push_back(glm::vec3{ radius * cos(theta + increment) * sin(phi + increment), radius * sin(theta + increment) * sin(phi + increment), radius * cos(phi + increment) });
			sphere.normals.push_back(glm::vec3{ cos(theta + increment) * sin(phi + increment), sin(theta + increment) * sin(phi + increment), cos(phi + increment) });
			sphere.texCoords.push_back(glm::vec2{1- (theta + increment) / (2 * M_PI), (phi + increment) / (1 * M_PI) });

			sphere.verts.push_back(glm::vec3{ radius * cos(theta + increment) * sin(phi), radius * sin(theta + increment) * sin(phi), radius * cos(phi) });
			sphere.normals.push_back(glm::vec3{ cos(theta + increment) * sin(phi), sin(theta + increment) * sin(phi), cos(phi) });
			sphere.texCoords.push_back(glm::vec2{1- (theta + increment) / (2 * M_PI), phi / (1 * M_PI) });
		}
	}
}

void applyShader(ShaderProgram& sp, glm::mat4& M, glm::mat4& P) {
	//glm::mat4 M = glm::mat4(1.0);
	//glm::mat4 V = camera.getView();
	//glm::mat4 P = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 1000.f);

	GLint location = glGetUniformLocation(sp, "light");
	glm::vec3 light = glm::vec3(0.0f, 0.0f, 1.0f);
	glUniform3fv(location, 1, glm::value_ptr(light));

	GLint uniMat = glGetUniformLocation(sp, "M");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(M));
	//uniMat = glGetUniformLocation(sp, "V");
	//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(V));
	uniMat = glGetUniformLocation(sp, "P");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(P));
}


int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(2160, 2160, "CPSC 453"); // can set callbacks at construction if desired


	GLDebug::enable();

	// CALLBACKS
	auto a4 = std::make_shared<Assignment4>();
	window.setCallbacks(a4);


	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.

	std::vector<glm::vec3> originQuad;
	originQuad.push_back(glm::vec3{-0.5, 0.5, 0}); // top-left
	originQuad.push_back(glm::vec3{-0.5, -0.5, 0}); // bottom-left
	originQuad.push_back(glm::vec3{0.5, 0.5, 0}); // top-right

	originQuad.push_back(glm::vec3{-0.5, -0.5, 0}); // bottom-left
	originQuad.push_back(glm::vec3{0.5, -0.5, 0}); // bottom-right
	originQuad.push_back(glm::vec3{0.5, 0.5, 0}); // top-right

	CPU_Geometry square;
	positiveZFace(originQuad, square);
	positiveXFace(originQuad, square);
	negativeZFace(originQuad, square);
	negativeXFace(originQuad, square);
	positiveYFace(originQuad, square);
	negativeYFace(originQuad, square);

	//square.cols.resize(square.verts.size(), glm::vec3{1.0, 0.0, 0.0});
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);
	colouredTriangles(square);


	/*for(auto i = square.verts.begin(); i < square.verts.end(); ++i) {
		std::cout << *i << std::endl;
	}*/

	GPU_Geometry quads;
	updateGPUGeometry(quads, square);

	GameObject Earth("images/2_no_clouds_16k.jpg", GL_LINEAR);
	GameObject Sun("images/8k_sun.jpg", GL_LINEAR);
	GameObject Moon("images/8k_moon.jpg", GL_LINEAR);
	GameObject Galaxy("images/8k_stars_milky_way.jpg", GL_LINEAR);


	//CPU_Geometry Earth;
	
	createSphere(0.1, 40, Earth.cgeom);
	createSphere(0.3, 40, Sun.cgeom);
	createSphere(0.3, 40, Moon.cgeom);
	createSphere(0.3, 40, Galaxy.cgeom);

	Earth.cgeom.cols.resize(Earth.cgeom.verts.size());
	Sun.cgeom.cols.resize(Sun.cgeom.verts.size());
	Moon.cgeom.cols.resize(Moon.cgeom.verts.size());
	Galaxy.cgeom.cols.resize(Galaxy.cgeom.verts.size());

	//Earth.transformationMatrix = glm::translate(Earth.transformationMatrix, glm::vec3(0.7f, 0.f, 0.f));
	//Earth.position += glm::vec3(0.7f, 0.f, 0.f);
	Earth.transformationMatrix = glm::scale(Earth.transformationMatrix, glm::vec3(0.1f));
	Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, glm::radians(90.f), glm::vec3(1.f, 0.f, 0.f));
	//reNormal(Earth.cgeom, glm::radians(90.f), 0);
	Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, glm::radians(-23.5f), glm::vec3(0.f, 1.f, 0.f));
	Moon.transformationMatrix = glm::rotate(Moon.transformationMatrix, glm::radians(6.7f), glm::vec3(0.f, 1.f, 0.f));
	glm::mat4 EarthDefault = Earth.transformationMatrix;

	Sun.transformationMatrix = glm::scale(Sun.transformationMatrix, glm::vec3(0.2f));
	glm::mat4 SunDefault = Sun.transformationMatrix;
	//Moon.transformationMatrix = glm::translate(Moon.transformationMatrix, glm::vec3(0.9f, 0.f, 0.f));
	Moon.transformationMatrix = glm::scale(Moon.transformationMatrix, glm::vec3(0.040f));
	glm::mat4 MoonDefault = Moon.transformationMatrix;
	Galaxy.transformationMatrix = glm::scale(Galaxy.transformationMatrix, glm::vec3(3.2f));

	//updateGPUGeometry(quads, Earth);
	updateGPUGeometry(Earth.ggeom, Earth.cgeom);
	updateGPUGeometry(Sun.ggeom, Sun.cgeom);
	updateGPUGeometry(Moon.ggeom, Moon.cgeom);
	updateGPUGeometry(Galaxy.ggeom, Galaxy.cgeom);

	glm::mat4 projection = glm::lookAt(glm::vec3(-0.2f, -0.2f, 0.60f),
		glm::vec3(0.0f, 0.0f, 0.0f),
		glm::vec3(0.0f, 1.0f, 0.0f));

	glm::mat4 temp = glm::mat4(1.0f);
	float angle = 0.f;
	float moonAngle = 0.f;
	float dx = 0;
	float dy = 0;
	float x = 0;
	float y = 0;
	float z = 0;
	float moonX = 0;
	float moonY = 0;
	float moonZ = 0;
	//Texture EarthTex("images/2_no_clouds_16k.jpg", GL_LINEAR);

	//Earth.cgeom.verts = glm::translate(glm::mat4(), glm::vec3(0.5f, 0, 0)) * Earth.verts;
	

	// RENDER LOOP
	while (!window.shouldClose()) {
		glfwPollEvents();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		shader.use();

		//a4->viewPipeline(shader);
		if (a4->buttonDown(GLFW_KEY_R)) {
			Earth.transformationMatrix = EarthDefault;
			Moon.transformationMatrix = MoonDefault;
			Sun.transformationMatrix = SunDefault;
			angle = 0; moonAngle = 0; x = 0; y = 0; z = 0; moonX = 0; moonY = 0; moonZ = 0;
		}
		if (!(a4->buttonDown(GLFW_KEY_SPACE))) {
			Sun.transformationMatrix = glm::rotate(Sun.transformationMatrix, -(float)(glfwGetTime() / 100000.0), glm::vec3(0.f, 1.f, 0.f));

			Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, -angle * 365.25f, glm::vec3(0.f, 0.f, 1.f));
			Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, glm::radians(23.5f), glm::vec3(0.f, 1.f, 0.f));
			Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, glm::radians(-90.f), glm::vec3(1.f, 0.f, 0.f));
			Earth.transformationMatrix = glm::scale(Earth.transformationMatrix, glm::vec3(1 / 0.1f));
			Earth.transformationMatrix = glm::translate(Earth.transformationMatrix, glm::vec3(-x, -z, -y));

			Moon.transformationMatrix = glm::rotate(Moon.transformationMatrix, -angle * 365.25f / 27, glm::vec3(0.f, 0.f, 1.f));
			Moon.transformationMatrix = glm::rotate(Moon.transformationMatrix, glm::radians(-6.7f), glm::vec3(0.f, 1.f, 0.f));
			Moon.transformationMatrix = glm::scale(Moon.transformationMatrix, glm::vec3(1 / 0.040f));
			Moon.transformationMatrix = glm::translate(Moon.transformationMatrix, glm::vec3(-x - moonX, -z - moonZ, -y - moonY));

			angle += 0.00001f;
			moonAngle += 27.0 * 0.00001f;
			x = 0.3 * cos(angle);
			y = 0.3 * sin(angle);
			z = 0.1 * sin(angle + M_PI / 4);
			moonX = 0.1 * cos(moonAngle);
			moonY = 0.1 * sin(moonAngle);
			moonZ = 0.03 * sin(moonAngle);


			Moon.transformationMatrix = glm::translate(Moon.transformationMatrix, glm::vec3(x + moonX, z + moonZ, y + moonY));
			Moon.transformationMatrix = glm::scale(Moon.transformationMatrix, glm::vec3(0.040f));
			Moon.transformationMatrix = glm::rotate(Moon.transformationMatrix, glm::radians(6.7f), glm::vec3(0.f, 1.f, 0.f));
			Moon.transformationMatrix = glm::rotate(Moon.transformationMatrix, angle * 365.25f / 27, glm::vec3(0.f, 0.f, 1.f));

			Earth.transformationMatrix = glm::translate(Earth.transformationMatrix, glm::vec3(x, z, y));
			Earth.transformationMatrix = glm::scale(Earth.transformationMatrix, glm::vec3(0.1f));
			Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, glm::radians(90.f), glm::vec3(1.f, 0.f, 0.f));
			Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, glm::radians(-23.5f), glm::vec3(0.f, 1.f, 0.f));
			Earth.transformationMatrix = glm::rotate(Earth.transformationMatrix, angle * 365.25f, glm::vec3(0.f, 0.f, 1.f));
		}


		//quads.bind();
		Earth.ggeom.bind();
		Sun.ggeom.bind();
		Moon.ggeom.bind();
		Galaxy.ggeom.bind();

		//glDrawArrays(GL_TRIANGLES, 0, GLsizei(square.verts.size()));
		applyShader(shader, Earth.transformationMatrix, projection);
		Earth.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(Earth.cgeom.verts.size()));
		Earth.texture.unbind();

		applyShader(shader, Sun.transformationMatrix, projection);
		Sun.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(Sun.cgeom.verts.size()));
		Sun.texture.unbind();

		applyShader(shader, Moon.transformationMatrix, projection);
		Moon.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(Moon.cgeom.verts.size()));
		Moon.texture.unbind();

		applyShader(shader, Galaxy.transformationMatrix, projection);
		Galaxy.texture.bind();
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(Galaxy.cgeom.verts.size()));
		Galaxy.texture.unbind();

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
