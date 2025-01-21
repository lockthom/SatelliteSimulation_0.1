#version 330 core

// Layout attributes
layout (location = 0) in vec3 aPos;   // Position: 1, 2, 3
layout (location = 1) in vec3 aColor; // Color: 4, 5, 6
layout (location = 2) in vec2 aTex;   // Texture: 7, 8

// Outputs
out vec3 color;
out vec2 texCoord;

// Inputs
uniform mat4 camMatrix;

void main() {

	gl_Position = camMatrix * vec4(aPos, 1.0);	

	// Outputs
	color = aColor;
	texCoord = aTex;
}
