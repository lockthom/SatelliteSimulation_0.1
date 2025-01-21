// General Imports
#include<iostream>        // Regular C++ stuff
#include<glad/glad.h>     // OpenGL preparation/ease of use
#include<GLFW/glfw3.h>    // OpenGL
#include<stb/stb_image.h> // Textures
#include<math.h>          // May not be strictly necessary

// openGL Matrix Stuff
#include<glm/glm.hpp>                  // Matrix Operations Overall (for now)
#include<glm/gtc/matrix_transform.hpp> // "
#include<glm/gtc/type_ptr.hpp>         // Allow pointers

// User defined classes
#include"Texture.h"     // Textures
#include"shaderClass.h" // Make shaders
#include"VAO.h"         // Vertex Array Object
#include"VBO.h"         // Vertex Buffer Object
#include"EBO.h"         // Index Buffer Object
#include"Camera.h"


// Window Parameters
const unsigned int windowWidth =  800;
const unsigned int windowHeight = 800;

int main() {

	// Initialize glfw
	glfwInit();

	// Clarify the version of openGL
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

	// What openGL profile do we want to use? (Packages of functions)
	// Core and Compatibility are the main options. Core the modern one.
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLfloat vertices[] = {
		// Pyramid
		-0.5f,  0.0f,  0.5f,     0.83f, 0.70f, 0.44f,    0.0f, 0.0f, // -x+z Base
		-0.5f,  0.0f, -0.5f,     0.83f, 0.70f, 0.44f,    5.0f, 0.0f, // -x-z Base
		 0.5f,  0.0f, -0.5f,     0.83f, 0.70f, 0.44f,    0.0f, 0.0f, // +x-z Base
		 0.5f,  0.0f,  0.5f,     0.83f, 0.70f, 0.44f,    5.0f, 0.0f, // +x+z Base
		 0.0f,  0.8f,  0.0f,     0.92f, 0.86f, 0.76f,    2.5f, 5.0f  // +y Tip
		// Square for textures (Verts, Colors, Tex. map)
		//-0.5f, -0.5f, 0.0f,     1.0f, 0.0f, 0.0f,    0.0f, 0.0f, // Lower left corner
		//-0.5f,  0.5f, 0.0f,     0.0f, 1.0f, 0.0f,    0.0f, 1.0f, // Upper left corner
		// 0.5f,  0.5f, 0.0f,     0.0f, 0.0f, 1.0f,    1.0f, 1.0f, // Upper right corner
		// 0.5f, -0.5f, 0.0f,     1.0f, 1.0f, 1.0f,    1.0f, 0.0f // Lower right corner
		// Zelda Triangle + Side Triangle ("XYZ")         ;/;/ RGBA Color
		//-0.5f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.8f, 0.3f,  0.02f, // Lower left corner
		// 0.5f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.8f, 0.3f,  0.02f, // Lower right corner
		// 0.0f,      0.5f * float(sqrt(3)) * 2 / 3, 0.0f,   1.0f, 0.6f,  0.32f, // Top corner
		//-0.5f / 2,  0.5f * float(sqrt(3)) / 6,     0.0f,   0.9f, 0.45f, 0.17f, // Middle left
		// 0.5f / 2,  0.5f * float(sqrt(3)) / 6,     0.0f,   0.9f, 0.45f, 0.17f, // Middle right
		// 0.0f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.8f, 0.3f,  0.02f, // Middle down
		// 0.75f,     0.5f * float(sqrt(3)) / 6,     0.0f,   0.0f, 1.0f,  0.0f, // Side triangle top
		// 1.0f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.0f, 0.0f,  1.0f  // Side triangle right
		//// Letter T
		//-0.1f, -0.4f, 0.0f, // Bottom left
		//-0.1f, 0.2f, 0.0f, // Mid left
		//-0.3f, 0.2f, 0.0f, // Overhand left
		//-0.3f, 0.4f, 0.0f, // Top left
		//0.3f, 0.4f, 0.0f, // Top right
		//0.3f, 0.2f, 0.0f, // Overhand right
		//0.1f, 0.2f, 0.0f, // Mid right
		//0.1f, -0.4f, 0.0f, // Bottom right
	};

	GLuint indices[] = {
		// Pyramid faces (Normal vector rotation outward)
		0, 1, 2, // 1/2 Base 
		0, 2, 3, // 2/2 Base
		0, 3, 4, // +z face
		2, 4, 3, // +x face
		1, 4, 2, // -z face
		0, 4, 1  // -x face
		//// Square for textures
		//0, 2, 1, // Upper triangle
		//0, 3, 2 // Lower triangle
		// Zelda + Side Indices
		//0, 3, 5, // Lower Left
		//3, 2, 4, // Lower Right
		//5, 4, 1, // Upper Triangle
		//1, 7, 6  // Side Triangle
		// Letter T indices
		//0, 1, 7, // Base left
		//7, 1, 6, // Base right
		//2, 3, 5, // Cap lower
		//5, 3, 4  // Cap Upper
	};

	// Now we generate a window.
	GLFWwindow* window = glfwCreateWindow(windowWidth, windowHeight, "Thomas Testing", NULL, NULL);
	
	// Make sure that if there are issues, glfw terminates.
	if (window == NULL) {

		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return EXIT_FAILURE;

	}

	// Not sure why this always fails. Glad is working?
	//if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	//{
	//	std::cout << "Failed to initialize GLAD" << std::endl;
	//	return EXIT_FAILURE;
	//}

	// Tell openGL we want to use the window.
	glfwMakeContextCurrent(window);

	// Load up glad.
	gladLoadGL();

	// Specify what openGL can "see" in the window
	glViewport(0, 0, windowWidth, windowHeight);

	// Declare shader with the vertex and fragment files
	Shader shaderProgram("default.vert", "default.frag");

	// Create a Vertex Array Object (VAO)
	VAO VAO1;
	VAO1.Bind();

	// Create Vertex Buffer Object (VBO) and Index Buffer (EBO)
	VBO VBO1(vertices, sizeof(vertices));
	EBO EBO1(indices, sizeof(indices));

	// Link info from VBO (position, color, texture) to VAO.
	// Layout value, number of values, type of values, stride to next set of values, offset from start.
	VAO1.LinkAttrib(VBO1, 0, 3, GL_FLOAT, 8 * sizeof(float), (void*)0);
	VAO1.LinkAttrib(VBO1, 1, 3, GL_FLOAT, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	VAO1.LinkAttrib(VBO1, 2, 2, GL_FLOAT, 8 * sizeof(float), (void*)(6 * sizeof(float)));

	// Unbind for consistency and to prevent unwanted edits
	VAO1.Unbind();
	VBO1.Unbind();
	EBO1.Unbind();

	// Texture
	// GL_RBG --> .jpg
	// GL_RBGA --> .png
	Texture knifeDuck("UW_Logo.png", GL_TEXTURE_2D, GL_TEXTURE0, GL_RGBA, GL_UNSIGNED_BYTE);
	knifeDuck.texUnit(shaderProgram, "tex0", 0);

	//// Pyramid rotation parameters
	//float rotation = 0.0f;
	//double prevTime = glfwGetTime();

	glEnable(GL_DEPTH_TEST);

	// Camera

	Camera camera(windowWidth, windowHeight, glm::vec3(0.0f, 0.0f, 2.0f));

	// Run until there is a reason to close it.
	while(!glfwWindowShouldClose(window)) {
	  
		//// Checking for new input
		//processInput(window);

		// Specify color of background
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);

		// Clear back buffer and depth buffer for new values.
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Start the shader program on the graphics card
		shaderProgram.Activate();

		camera.Inputs(window);
		camera.Matrix(45.0f, 0.1f, 100.0f, shaderProgram, "camMatrix");

		/*double crntTime = glfwGetTime();
		if (crntTime - prevTime >= 1 / 60) {

			rotation += 2.0f;
			prevTime = crntTime;

		}*/
		
		// Bind texture
		knifeDuck.Bind();

		// Bind vertex array
		VAO1.Bind();

		// Draw! (12 for triangle, 6 for texture square)
		glDrawElements(GL_TRIANGLES, sizeof(indices)/sizeof(int), GL_UNSIGNED_INT, 0);

		// Swap buffers
		glfwSwapBuffers(window);

		// "Take care of all glfw events"
		glfwPollEvents();

	}


	// Clean up user info
	VAO1.Delete();
	VBO1.Delete();
	EBO1.Delete();
	knifeDuck.Delete();
	shaderProgram.Delete();

	// Clean up generally
	glfwDestroyWindow(window);
	glfwTerminate();

	return EXIT_SUCCESS;

}