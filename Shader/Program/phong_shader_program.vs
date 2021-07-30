#version 330 core
layout (location = 0) in vec3 a_position;
layout (location = 1) in vec2 a_texcoord;

out vec3 fragPos;
out vec2 texcoord;
out vec4 lightSpacePos;

uniform mat4 projectMatrix;
uniform mat4 viewMatrix;
uniform mat4 lightSpaceMatrix;

void main()
{
	gl_Position = projectMatrix * viewMatrix * vec4(a_position, 1.0);
	fragPos = a_position;
	texcoord = a_texcoord;
	lightSpacePos = lightSpaceMatrix * vec4(a_position, 1.0);
} 