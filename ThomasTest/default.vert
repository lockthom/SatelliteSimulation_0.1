#version 330 core

// Position
layout (location = 0) in vec3 aPos;

// Color
layout (location = 1) in vec3 aColor;

// Provide color for the fragment shader
out vec3 color;

// Scale the vertices.
uniform float scale;

// World Matrices
uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main() {

	// Scale and then apply all world matrices, in order.
	gl_Position = proj * view * model * vec4(aPos.x + aPos.x * scale, aPos.y + aPos.y * scale, aPos.z + aPos.z * scale, 1.0);
	
	color = aColor;

}
