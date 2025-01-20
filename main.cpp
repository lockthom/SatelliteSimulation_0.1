// General Imports
#include<iostream>
#include<glad/glad.h>
#include<GLFW/glfw3.h>
#include<math.h> // May notbe strictly necessary

// openGL Matrix Stuff
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>

// User defined classes
#include"shaderClass.h"
#include"VAO.h"
#include"VBO.h"
#include"EBO.h"

const unsigned int windowWidth =  800;
const unsigned int windowHeight = 800;


void processInput(GLFWwindow* window){

	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

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
		// Zelda Triangle + Side Triangle ("XYZ")         ;/;/ RGBA Color
		-0.5f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.8f, 0.3f,  0.02f, // Lower left corner
		 0.5f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.8f, 0.3f,  0.02f, // Lower right corner
		 0.0f,      0.5f * float(sqrt(3)) * 2 / 3, 0.0f,   1.0f, 0.6f,  0.32f, // Top corner
		-0.5f / 2,  0.5f * float(sqrt(3)) / 6,     0.0f,   0.9f, 0.45f, 0.17f, // Middle left
		 0.5f / 2,  0.5f * float(sqrt(3)) / 6,     0.0f,   0.9f, 0.45f, 0.17f, // Middle right
		 0.0f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.8f, 0.3f,  0.02f, // Middle down
		 0.75f,     0.5f * float(sqrt(3)) / 6,     0.0f,   0.0f, 1.0f,  0.0f, // Side triangle top
		 1.0f,     -0.5f * float(sqrt(3)) / 3,     0.0f,   0.0f, 0.0f,  1.0f  // Side triangle right
		// Letter T
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
		// Zelda + Side Indices
		0, 3, 5, // Lower Left
		3, 2, 4, // Lower Right
		5, 4, 1, // Upper Triangle
		1, 7, 6  // Side Triangle
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

	// Declare shader
	Shader shaderProgram("default.vert", "default.frag");

	// Create a Vertex Array Object (VAO)
	VAO VAO1;
	VAO1.Bind();

	// Create Vertex Buffer Object (VBO) and Index Buffer (EBO)
	VBO VBO1(vertices, sizeof(vertices));
	EBO EBO1(indices, sizeof(indices));

	// Link info from VBO (position and color) to VAO.
	VAO1.LinkAttrib(VBO1, 0, 3, GL_FLOAT, 6 * sizeof(float), (void*)0);
	VAO1.LinkAttrib(VBO1, 1, 3, GL_FLOAT, 6 * sizeof(float), (void*)(3 * sizeof(float)));

	// Unbind for consistency
	VAO1.Unbind();
	VBO1.Unbind();
	EBO1.Unbind();

	// Create ID for the "uniform" (a value) in the shader for scaling.
	GLuint uniID = glGetUniformLocation(shaderProgram.ID, "scale");

	// Bring buffer with color to display
	glfwSwapBuffers(window);

	// Run until there is a reason to close it.
	while(!glfwWindowShouldClose(window)) {
	
		// Checking for new input
		processInput(window);

		// Specify color of background
		glClearColor(0.07f, 0.13f, 0.17f, 1.0f);

		// Clear back buffer, assign new color
		glClear(GL_COLOR_BUFFER_BIT);

		// Start shader, point to scaling value.
		shaderProgram.Activate();
		glUniform1f(uniID, 0.0f);

		// Initialize relevant world matrices
		glm::mat4 model = glm::mat4(1.0f);
		glm::mat4 view = glm::mat4(1.0f);
		glm::mat4 proj = glm::mat4(1.0f);

		// Fill relevant world matrices
		//                                      X     Y      Z
		view = glm::translate(view, glm::vec3(0.0f, 0.0f, -2.0f));
		proj = glm::perspective(glm::radians(45.0f), (float)(windowWidth / windowHeight), 0.1f, 100.0f);

		// Connect world matrices to shader
		int modelLoc = glGetUniformLocation(shaderProgram.ID, "model");
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

		int viewLoc = glGetUniformLocation(shaderProgram.ID, "view");
		glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));

		int projLoc = glGetUniformLocation(shaderProgram.ID, "proj");
		glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));
		
		// Bind vertex array
		VAO1.Bind();

		// Draw!
		glDrawElements(GL_TRIANGLES, 12, GL_UNSIGNED_INT, 0);

		// Swap buffers
		glfwSwapBuffers(window);

		// "Take care of all glfw events"
		glfwPollEvents();

	}


	// Clean up user info
	VAO1.Delete();
	VBO1.Delete();
	EBO1.Delete();
	shaderProgram.Delete();

	// Clean up generally
	glfwDestroyWindow(window);
	glfwTerminate();

	return EXIT_SUCCESS;

}