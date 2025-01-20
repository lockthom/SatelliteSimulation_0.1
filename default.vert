#version 330 core

// Position format in the buffer (first 3 values)
layout (location = 0) in vec3 aPos;

// Color format in the buffer (second 3 values, or 4, 5, and 6)
layout (location = 1) in vec3 aColor;

// Output color. Goes to the fragment shader
out vec3 color;

// Scale the vertices.
// This is basically a constant we must pass into GLSL and must use.
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
