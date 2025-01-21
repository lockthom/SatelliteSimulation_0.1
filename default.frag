
// Specify version
#version 330 core

// Output of the fragment shader is a RGBA vector
out vec4 FragColor;

// This input to the fragment shader is a RGB vector (from vertex shader)
in vec3 color;

// This input includes the texture values (from the vertex shader)
in vec2 texCoord;

uniform sampler2D tex0;

void main()
{
	//FragColor = vec4(color, 1.0f);
	FragColor = texture(tex0, texCoord);
}