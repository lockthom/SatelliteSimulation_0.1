
// Specify version
#version 330 core

// Output of the fragment shader is a RGBA vector
out vec4 FragColor;

// Input to the fragment shader is a RGB vector
in vec3 color;

void main()
{
	FragColor = vec4(color, 1.0f);
}