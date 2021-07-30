#version 330 core
layout (location = 0) in vec3 a_position;

uniform mat4 lightSpaceMatrix;
out vec4 lightSpacePos;

void main()
{
    lightSpacePos = lightSpaceMatrix * vec4(a_position, 1.0);
    gl_Position = lightSpacePos;
}