#include"Camera.h"

Camera::Camera(int width, int height, glm::vec3 position) {

	Camera::width = width;
	Camera::height = height;
	Position = position;
}

void Camera::Matrix(float FOVdeg, float nearPlane, float farPlane, Shader& shader, const char* uniform) {

	glm::mat4 view = glm::mat4(1.0f);
	glm::mat4 projection = glm::mat4(1.0f);

	view = glm::lookAt(Position, Position + Orientation, Up);
	projection = glm::perspective(glm::radians(FOVdeg), (float)(height / width), nearPlane, farPlane);

	glUniformMatrix4fv(glGetUniformLocation(shader.ID, uniform), 1, GL_FALSE, glm::value_ptr(projection * view));

}

void Camera::Inputs(GLFWwindow* window) {

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {

		// Forward motion
		Position += speed * Orientation;

	}
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {

		// Move to camera left
		Position += speed * -glm::normalize(glm::cross(Orientation, Up));

	}
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {

		// Backwards motion
		Position += speed * -Orientation;

	}
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {

		// Move to camera right
		Position += speed * glm::normalize(glm::cross(Orientation, Up));

	}
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {

		// Vertical upwards motion
		Position += speed * Up;

	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) {

		// Vertical downwards motion
		Position += speed * -Up;

	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {

		// Increase speed while shift is pressed
		speed = 0.04f;

	}
	else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_RELEASE) {

		// Return to normal speed upon release of shift
		speed = 0.01f;

	}
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {

		// Close window when escape is pressed
		glfwSetWindowShouldClose(window, true);

	}

	// Mouse inputs

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {

		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

		if (firstClick) {

			glfwSetCursorPos(window, (width / 2), (height / 2));
			firstClick = false;
		}

		double mouseX;
		double mouseY;
		glfwGetCursorPos(window, &mouseX, &mouseY);

		float rotX = sensitivity * (float)(mouseY - (height / 2)) / height;
		float rotY = sensitivity * (float)(mouseX - (height / 2)) / height;

		glm::vec3 newOrientation = glm::rotate(Orientation, glm::radians(-rotX), glm::normalize(glm::cross(Orientation, Up)));

		if (!((glm::angle(newOrientation, Up) <= glm::radians(5.0f)) or (glm::angle(newOrientation, -Up) <= glm::radians(5.0f)))) {

			Orientation = newOrientation;

		}

		Orientation = glm::rotate(Orientation, glm::radians(-rotY), Up);

		glfwSetCursorPos(window, (width / 2), (height / 2));

	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) {

		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
		firstClick = true;

	}
}