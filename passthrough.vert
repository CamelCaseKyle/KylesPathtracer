#version 330 core
precision highp int;
precision highp float;

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec2 position;
layout(location = 1) in vec2 texcoord;

// Output data; will be interpolated for each fragment.
smooth out vec2 uv;

uniform sampler2D iChannel0;
uniform sampler2D iChannel1; 
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform vec2 iResolution;
uniform int iFrame;
uniform float iTime;
uniform vec3 loc;
uniform vec3 vel;
uniform vec3 iMouse;
uniform vec3 orient;

void main() {
	uv = texcoord;
	gl_Position = vec4(position, 0., 1.);
}
