#version 330 core

layout(location = 0) out vec4 fragColor;

in vec2 uv;

uniform sampler2D iChannel0;
uniform sampler2D iChannel1; 
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform vec2 iResolution;
uniform int iFrame;
uniform float iGlobalTime;
uniform vec3 loc;
uniform vec3 vel;
uniform vec3 iMouse;
uniform vec3 orient;

void main() {
	vec4 col = texture2D(iChannel0, uv * .5 + .5);
	fragColor = col / col.a;
}